/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/

#include "dfChemistryModel.H"
#include "CanteraMixture.H"
#include "hePsiThermo.H"

#ifdef USE_PYTORCH
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h> //used to convert
#endif

#ifdef USE_LIBTORCH
#include <torch/script.h>
#include "DNNInferencer.H"
#endif 

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "PstreamGlobals.H"
#include "CombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#ifdef USE_PYTORCH
    pybind11::scoped_interpreter guard{};//start python interpreter
#endif
    #define NO_CONTROL

    #include "postProcess.H"

    // #include "setRootCaseLists.H"
    #include "listOptions.H"
    #include "setRootCase2.H"
    #include "listOutput.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    double time_monitor_flow=0;
    double time_monitor_chem=0;
    double time_monitor_Y=0;
    double time_monitor_AMR=0;
    clock_t start, end;

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    dimensionedScalar refCri("refCri", dimensionSet(1, -4, 0, 0, 0), 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        //used for AMR
        refCri = max(mag(fvc::grad(rho)));
        tmp<volScalarField> tmagGradrho = mag(fvc::grad(rho));
        volScalarField normalisedGradrho
        (
            "normalisedGradrho",
            tmagGradrho()/refCri
        );
        normalisedGradrho.writeOpt() = IOobject::AUTO_WRITE;
        tmagGradrho.clear();

        if (!LTS)
        {
            #include "setDeltaT.H"
            runTime++;

            // Do any mesh changes
            start = std::clock();
            mesh.update();
            end = std::clock();
            time_monitor_AMR += double(end - start) / double(CLOCKS_PER_SEC);
            
        }

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            phiv_pos -= mesh.phi();
            phiv_neg -= mesh.phi();
        }

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
            runTime++;
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
        }

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        #include "rhoEqn.H"

        start = std::clock();
        // --- Solve momentum
        #include "rhoUEqn.H"
        end = std::clock();
        time_monitor_flow += double(end - start) / double(CLOCKS_PER_SEC);


        // --- Solve species
        #include "rhoYEqn.H"

        // --- Solve energy
       #include "rhoEEqn.H"

        turbulence->correct();

        runTime.write();

        Info<< "MonitorTime_chem = " << time_monitor_chem << " s" << nl << endl;
        Info<< "MonitorTime_Y = " << time_monitor_Y << " s" << nl << endl;
        Info<< "MonitorTime_flow = " << time_monitor_flow << " s" << nl << endl;
        Info<< "MonitorTime_AMR = " << time_monitor_AMR << " s" << nl << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
