/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "UList.H"

#include "hexRefData.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "dfRefinementHistory.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexRefData::hexRefData(const IOobject& io)
{
    {
        IOobject rio(io);
        rio.rename("cellLevel");
        // bool haveFile = returnReduce(rio.headerOk(), orOp<bool>());
        bool haveFile = returnReduce(rio.typeHeaderOk<labelIOList>(), orOp<bool>());
        if (haveFile)
        {
            Info<< "Reading hexRef data : " << rio.name() << endl;
            cellLevelPtr_.reset(new labelIOList(rio));
        }
    }
    {
        IOobject rio(io);
        rio.rename("pointLevel");
        // bool haveFile = returnReduce(rio.headerOk(), orOp<bool>());
        bool haveFile = returnReduce(rio.typeHeaderOk<labelIOList>(), orOp<bool>());
        if (haveFile)
        {
            Info<< "Reading hexRef data : " << rio.name() << endl;
            pointLevelPtr_.reset(new labelIOList(rio));
        }
    }
    {
        IOobject rio(io);
        rio.rename("level0Edge");
        // bool haveFile = returnReduce(rio.headerOk(), orOp<bool>());
        bool haveFile = returnReduce(rio.typeHeaderOk<labelIOList>(), orOp<bool>());
        if (haveFile)
        {
            Info<< "Reading hexRef data : " << rio.name() << endl;
            level0EdgePtr_.reset(new uniformDimensionedScalarField(rio));
        }
    }
    {
        IOobject rio(io);
        rio.rename("dfRefinementHistory");
        // bool haveFile = returnReduce(rio.headerOk(), orOp<bool>());
        bool haveFile = returnReduce(rio.typeHeaderOk<labelIOList>(), orOp<bool>());
        if (haveFile)
        {
            Info<< "Reading hexRef data : " << rio.name() << endl;
            refHistoryPtr_.reset(new dfRefinementHistory(rio));
        }
    }
}


Foam::hexRefData::hexRefData
(
    const IOobject& io,
    const hexRefData& data,
    const labelList& cellMap,
    const labelList& pointMap
)
{
    if (data.cellLevelPtr_.valid())
    {
        IOobject rio(io);
        rio.rename(data.cellLevelPtr_().name());

        cellLevelPtr_.reset
        (
            new labelIOList
            (
                rio,
                UIndirectList<label>(data.cellLevelPtr_(), cellMap)()
            )
        );
    }
    if (data.pointLevelPtr_.valid())
    {
        IOobject rio(io);
        rio.rename(data.pointLevelPtr_().name());

        pointLevelPtr_.reset
        (
            new labelIOList
            (
                rio,
                UIndirectList<label>(data.pointLevelPtr_(), pointMap)()
            )
        );
    }
    if (data.level0EdgePtr_.valid())
    {
        IOobject rio(io);
        rio.rename(data.level0EdgePtr_().name());

        level0EdgePtr_.reset
        (
            new uniformDimensionedScalarField(rio, data.level0EdgePtr_())
        );
    }
    if (data.refHistoryPtr_.valid())
    {
        IOobject rio(io);
        rio.rename(data.refHistoryPtr_().name());

        refHistoryPtr_ = data.refHistoryPtr_().clone(rio, cellMap);
    }
}


Foam::hexRefData::hexRefData
(
    const IOobject& io,
    const UPtrList<const labelList>& cellMaps,
    const UPtrList<const labelList>& pointMaps,
    const UPtrList<const hexRefData>& procDatas
)
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(io.db());

    // cellLevel

    if (procDatas[0].cellLevelPtr_.valid())
    {
        IOobject rio(io);
        rio.rename(procDatas[0].cellLevelPtr_().name());

        cellLevelPtr_.reset(new labelIOList(rio, mesh.nCells()));
        labelList& cellLevel = cellLevelPtr_();

        forAll(procDatas, procI)
        {
            const labelList& procCellLevel = procDatas[procI].cellLevelPtr_();
            UIndirectList<label>(cellLevel, cellMaps[procI]) = procCellLevel;
        }
    }


    // pointLevel

    if (procDatas[0].pointLevelPtr_.valid())
    {
        IOobject rio(io);
        rio.rename(procDatas[0].pointLevelPtr_().name());

        pointLevelPtr_.reset(new labelIOList(rio, mesh.nPoints()));
        labelList& pointLevel = pointLevelPtr_();

        forAll(procDatas, procI)
        {
            const labelList& procPointLevel = procDatas[procI].pointLevelPtr_();
            UIndirectList<label>(pointLevel, pointMaps[procI]) = procPointLevel;
        }
    }


    // level0Edge

    if (procDatas[0].level0EdgePtr_.valid())
    {
        IOobject rio(io);
        rio.rename(procDatas[0].level0EdgePtr_().name());

        level0EdgePtr_.reset
        (
            new uniformDimensionedScalarField
            (
                rio,
                procDatas[0].level0EdgePtr_()
            )
        );
    }


    // refinementHistory

    if (procDatas[0].refHistoryPtr_.valid())
    {
        IOobject rio(io);
        rio.rename(procDatas[0].refHistoryPtr_().name());

        UPtrList<const dfRefinementHistory> procRefs(procDatas.size());
        forAll(procDatas, i)
        {
            procRefs.set(i, &procDatas[i].refHistoryPtr_());
        }

        refHistoryPtr_.reset
        (
            new dfRefinementHistory
            (
                rio,
                cellMaps,
                procRefs
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hexRefData::~hexRefData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hexRefData::sync(const IOobject& io)
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(io.db());

    bool hasCellLevel = returnReduce(cellLevelPtr_.valid(), orOp<bool>());
    if (hasCellLevel && !cellLevelPtr_.valid())
    {
        IOobject rio(io);
        rio.rename("cellLevel");
        rio.readOpt() = IOobject::NO_READ;
        cellLevelPtr_.reset(new labelIOList(rio, labelList(mesh.nCells(), 0)));
    }

    bool hasPointLevel = returnReduce(pointLevelPtr_.valid(), orOp<bool>());
    if (hasPointLevel && !pointLevelPtr_.valid())
    {
        IOobject rio(io);
        rio.rename("pointLevel");
        rio.readOpt() = IOobject::NO_READ;
        pointLevelPtr_.reset
        (
            new labelIOList(rio, labelList(mesh.nPoints(), 0))
        );
    }

    bool hasLevel0Edge = returnReduce(level0EdgePtr_.valid(), orOp<bool>());
    if (hasLevel0Edge)
    {
        // Get master length
        scalar masterLen = level0EdgePtr_().value();
        Pstream::scatter(masterLen);
        if (!level0EdgePtr_.valid())
        {
            IOobject rio(io);
            rio.rename("level0Edge");
            rio.readOpt() = IOobject::NO_READ;
            level0EdgePtr_.reset
            (
                new uniformDimensionedScalarField
                (
                    rio,
                    dimensionedScalar("zero", dimLength, masterLen)
                )
            );
        }
    }

    bool hasHistory = returnReduce(refHistoryPtr_.valid(), orOp<bool>());
    if (hasHistory && !refHistoryPtr_.valid())
    {
        IOobject rio(io);
        rio.rename("dfRefinementHistory");
        rio.readOpt() = IOobject::NO_READ;
        refHistoryPtr_.reset(new dfRefinementHistory(rio, mesh.nCells(), true));
    }
}


void Foam::hexRefData::distribute(const mapDistributePolyMesh& map)
{
    if (cellLevelPtr_.valid())
    {
        map.cellMap().distribute(cellLevelPtr_());
    }
    if (pointLevelPtr_.valid())
    {
        map.pointMap().distribute(pointLevelPtr_());
    }

    // No need to distribute the level0Edge

    if (refHistoryPtr_.valid() && refHistoryPtr_().active())
    {
        refHistoryPtr_().distribute(map);
    }
}


bool Foam::hexRefData::write() const
{
    bool ok = true;
    if (cellLevelPtr_.valid())
    {
        ok = ok && cellLevelPtr_().write();
    }
    if (pointLevelPtr_.valid())
    {
        ok = ok && pointLevelPtr_().write();
    }
    if (level0EdgePtr_.valid())
    {
        ok = ok && level0EdgePtr_().write();
    }
    if (refHistoryPtr_.valid())
    {
        ok = ok && refHistoryPtr_().write();
    }
    return ok;
}


// ************************************************************************* //
