surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(turbulence->muEff())*mesh.magSf()*fvc::snGrad(U)
      + (mesh.Sf() & fvc::interpolate(tauMC))
    )
    & (a_pos*U_pos + a_neg*U_neg)
);

solve
(
      fvm::ddt(rhoE)
    + fvc::div(phiEp)
    - fvc::div(sigmaDotU)
);

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();

ha = e + p/rho;
chemistry->correctThermo(); // before this, we must update ha = e + p/rho

rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

if (!inviscid)
{
    fvScalarMatrix eEqn
    (
           fvm::ddt(rho, e) - fvc::ddt(rho, e)
        // alpha in deepflame is considered to calculate h by default (kappa/Cp), so multiply gamma to correct alpha
         - fvm::laplacian(turbulence->alphaEff()*thermo.gamma(), e)
         + diffAlphaD
      ==
           fvc::div(hDiffCorrFlux)
    );

    eEqn.solve("e");

    ha = e + p/rho;
    chemistry->correctThermo();
    rhoE = rho*(e + 0.5*magSqr(U));
}

Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

p.ref() = rho()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
