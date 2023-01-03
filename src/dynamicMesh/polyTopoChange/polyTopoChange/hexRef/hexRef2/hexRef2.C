/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "hexRef2.H"

#include "emptyPolyPatch.H"

#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef2, 0);
    addToRunTimeSelectionTable(hexRef, hexRef2, mesh);
    addToRunTimeSelectionTable(hexRef, hexRef2, levelsHist);
    addToRunTimeSelectionTable(hexRef, hexRef2, levels);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Check whether pointi is an anchor on celli.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRef2::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label celli,
    const label facei,
    const label pointi
) const
{
    return 1;
}

Foam::label Foam::hexRef2::myGetAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& edgeMidPoint,
    const label celli,
    const label facei,
    const label pointi
) const
{
    
    if (cellAnchorPoints[celli].size())
    {
        
        label index = 0;

        const face& f = mesh_.faces()[facei];
        label nextPointi = -1;

        forAll(f,fp)
        {
            nextPointi = f[fp];
            if( pointi == f[f.fcIndex(fp)] || pointi == f[f.rcIndex(fp)] )
            {
                label edgeJ = meshTools::findEdge (mesh_, pointi, nextPointi);
                if(edgeMidPoint[edgeJ] >= 0)
                {
                    // const edge& ei = mesh_.edges()[edgeJ];
                    // Pout<< "find anchor cell edge is :" << ei<< endl;
                    break;
                }
            }
        }
        if(pointi < nextPointi)
        {
            index = 0;
        }
        else
        {
            index = 1;
        }

        //use for check
        // Pout<< "go this way : 1" << endl;
        // Pout<< "pointi is:" << pointi << endl;
        // Pout<< "nextPointi is:" << nextPointi << endl;
        // Pout<< "facei is: "<<f << endl;
        // Pout<< "index is:" << index << endl;

        return cellAddedCells[celli][index];

        // Problem.
        dumpCell(celli);
        Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
            << endl;

        FatalErrorInFunction
            << "Could not find point " << pointi
            << " in the anchorPoints for cell " << celli << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        //use for check
        // Pout<< "go this way : 3" << endl;
        return celli;
    }
}

// not be used in hexRef2
Foam::label Foam::hexRef2::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& edgeMidPoint,
    const label celli,
    const label facei,
    const bool faceOrder,
    const label edgeMidPointi,
    const label anchorPointi,
    const label faceMidPointi,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    return -1;
}


//not be used in hexRef2
void Foam::hexRef2::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label celli,

    polyTopoChange& meshMod
) const
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef2::hexRef2(const polyMesh& mesh, const bool readHistory)
:
    hexRef(mesh, readHistory)
{}


// Construct from components
Foam::hexRef2::hexRef2
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const dfRefinementHistory& history,
    const scalar level0Edge
)
:
    hexRef(mesh, cellLevel, pointLevel, history, level0Edge)
{}


// Construct from components
Foam::hexRef2::hexRef2
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    hexRef(mesh, cellLevel, pointLevel, level0Edge)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Top level driver to insert topo changes to do all refinement.
Foam::labelListList Foam::hexRef2::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        //checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();


    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, celli)
    {
        newCellLevel.append(cellLevel_[celli]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointi)
    {
        newPointLevel.append(pointLevel_[pointi]);
    }


    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Allocating " << cellLabels.size() << " cell midpoints."
            << endl;
    }


    // Mid point per refined cell.
    // -1 : not refined
    // >=0: label of mid point.
    // just used to sign the cell that need to be split
    labelList cellMidPoint(mesh_.nCells(), -1);

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];
        
        cellMidPoint[celli] = 12345;
    }

    boolList isDivisibleFace(mesh_.nFaces(), false);
    boolList isDivisibleEdge(mesh_.nEdges(), false);

    // the last face and first face cannot cannot be divised
    label firstBFace = -1;
    label lastBFace = -1;
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const label & patchID = mesh_.boundaryMesh().whichPatch(facei);

        if (mesh_.boundaryMesh()[patchID].type() == "empty")
        {
            isDivisibleFace[facei] = true;
            const labelList& fEdges = mesh_.faceEdges(facei);

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];
                isDivisibleEdge[edgeI] = true;
            }
        }
        else
        {
            if( firstBFace == -1 )
            {
                firstBFace = facei;
            }
            else
            {
                lastBFace = facei;
            }
        }
    }

    isDivisibleFace[firstBFace] = false;
    const labelList& fbfEdges = mesh_.faceEdges(firstBFace);
    forAll(fbfEdges, i)
    {
        label edgeI = fbfEdges[i];
        isDivisibleEdge[edgeI] = false;
    }
    
    isDivisibleFace[lastBFace] = false;
    const labelList& lbfEdges = mesh_.faceEdges(lastBFace);
    forAll(lbfEdges, i)
    {
        label edgeI = lbfEdges[i];
        isDivisibleEdge[edgeI] = false;
    }

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        isDivisibleFace[facei] = false;
        const labelList& fEdges = mesh_.faceEdges(facei);

        forAll(fEdges, i)
        {
                label edgeI = fEdges[i];
                isDivisibleEdge[edgeI] = false;
        }
    }

    // use for check
    // forAll(cellLabels, i)
    // {
    //     label celli = cellLabels[i];
    //     if( cellMidPoint[celli] >= 0)
    //     {
    //         // label nface = 0;
    //         label nedges = 0;
    //         const labelList& cEdges = mesh_.cellEdges(celli);
    //         forAll(cEdges, i)
    //         {
    //             label edgeI = cEdges[i];
    //             if(isDivisibleEdge[edgeI])
    //             {
    //                 nedges++;
    //             }
    //         }
    //         Pout<< "each cell edges need to split:"
    //         << nedges << endl;
            
    //     }
    // }

    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellLabels.size());

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                splitCells.insert(celli);
            }
        }

        Pout<< "hexRef2::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }



    // Split edges
    // ~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over cells to be refined or edges?
    forAll(cellMidPoint, celli)
    {
        // label n = 0;  //use to check
        if (cellMidPoint[celli] >= 0)
        {
            const labelList& cEdges = mesh_.cellEdges(celli);

            // use to check
            // Pout<< "edges' number for a cell need to split : "<< cEdges.size()<<endl;
            // Pout << "edges are :" << endl; 
                        
            // forAll(cEdges, i)
            // {
            //     label edgeI = cEdges[i];
            //     const edge& e = mesh_.edges()[edgeI];
            //     Pout << e << endl; 
            // }
            

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];
                const edge& e = mesh_.edges()[edgeI];
                if
                (
                    isDivisibleEdge[edgeI]
                 && pointLevel_[e[0]] <= cellLevel_[celli]
                 && pointLevel_[e[1]] <= cellLevel_[celli]
                )
                {
                    edgeMidPoint[edgeI] = 12345;    // mark need for splitting
                    //n++;  //use to check
                }
            }
        }

        //use to check
        // if(n > 0)
        // {
        //     Pout<< "edges need to split:"
        //     << n << endl;
        // }
        
    }

    // Synchronize edgeMidPoint across coupled patches. Take max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
        );


        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    )
                  + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexRef2::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Calculate face level
    // ~~~~~~~~~~~~~~~~~~~~
    // (after splitting)

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label facei = 0; facei < mesh_.nFaces(); facei++)
    {
        faceAnchorLevel[facei] = faceLevel(facei);
    }

    // -1  : no need to split face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);


    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    // for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    // {
    //     if (faceAnchorLevel[facei] >= 0)
    //     {
    //         label own = mesh_.faceOwner()[facei];
    //         label ownLevel = cellLevel_[own];
    //         label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

    //         label nei = mesh_.faceNeighbour()[facei];
    //         label neiLevel = cellLevel_[nei];
    //         label newNeiLevel = neiLevel + (cellMidPoint[nei] >= 0 ? 1 : 0);

    //         if(isDivisibleFace[facei])
    //         {
    //             // if
    //             // (
    //             //     newOwnLevel > faceAnchorLevel[facei]
    //             // || newNeiLevel > faceAnchorLevel[facei]
    //             // )
    //             // {
    //                 faceMidPoint[facei] = 12345;    // mark to be split
    //             // }
    //         }
    //     }
    // }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceAnchorLevel[facei] >= 0)
            {
                //useless in hexRef2
                // label own = mesh_.faceOwner()[facei];
                // label ownLevel = cellLevel_[own];
                // label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

                if(isDivisibleFace[facei])
                {
                    faceMidPoint[facei] = 12345;    // mark to be split
                }
                
            }
        }
    }

    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );

    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~
    // no need to introduce face mid points

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", cellLabels.size());

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                splitFaces.insert(facei);
            }
        }

        Pout<< "hexRef2::setRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - cellMidPoint >= 0 : cell needs to be split
    // - faceMidPoint >= 0 : face needs to be split
    // - edgeMidPoint >= 0 : edge needs to be split
    // - isDivisibleFace true :face needs to be split in two



    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());

    {
        labelList nAnchorPoints(mesh_.nCells(), 0);

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                cellAnchorPoints[celli].setSize(8);
            }
        }

        forAll(pointLevel_, pointi)
        {
            const labelList& pCells = mesh_.pointCells(pointi);
           
            forAll(pCells, pCelli)
            {
                label celli = pCells[pCelli];

                if
                (
                    cellMidPoint[celli] >= 0
                 && pointLevel_[pointi] <= cellLevel_[celli]
                )
                {
                    if (nAnchorPoints[celli] == 8)
                    {
                        dumpCell(celli);
                        FatalErrorInFunction
                            << "cell " << celli
                            << " of level " << cellLevel_[celli]
                            << " uses more than 8 points of equal or"
                            << " lower level" << nl
                            << "Points so far:" << cellAnchorPoints[celli]
                            << abort(FatalError);
                    }

                    cellAnchorPoints[celli][nAnchorPoints[celli]++] = pointi;
                }
            }
        }

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                if (nAnchorPoints[celli] != 8)
                {
                    dumpCell(celli);

                    const labelList& cPoints = mesh_.cellPoints(celli);

                    FatalErrorInFunction
                        << "cell " << celli
                        << " of level " << cellLevel_[celli]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << UIndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Adding cells (1 per anchorPoint)"
            << endl;
    }

    // Per cell the 1 added cells (+ original cell)
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellAnchorPoints, celli)
    {
        const labelList& cAnchors = cellAnchorPoints[celli];

        if (cAnchors.size() == 8)
        {
            labelList& cAdded = cellAddedCells[celli];
            cAdded.setSize(2);

            // Original cell at 0
            cAdded[0] = celli;

            // Update cell level
            newCellLevel[celli] = cellLevel_[celli]+1;


            for (label i = 1; i < 2; i++)
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                                 // master point
                        -1,                                 // master edge
                        -1,                                 // master face
                        celli,                              // master cell
                        mesh_.cellZones().whichZone(celli)  // zone for cell
                    )
                );

                newCellLevel(cAdded[i]) = cellLevel_[celli]+1;
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but get new owner/neighbour
    // 3. new internal faces inside split cells.

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, i)
                {
                    label facei = cFaces[i];
                    if(!isDivisibleFace[facei])
                    {
                        affectedFace.set(cFaces[i]);
                    } 
                }
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef2::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, facei)
    {
        if (faceMidPoint[facei] >= 0)
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[facei];

            // Has original facei been used (one face added, original gets
            // modified)
            bool modifiedFace = false;

            //useless in hexRef2
            // label anchorLevel = faceAnchorLevel[facei];

            if(isDivisibleFace[facei])
            {
                face newFace(2);

                // maybe need to adjust the point order
                // DynamicList<label> faceVerts(4);              

                forAll(f, fp)
                {
                    //maybe need to put inside (but nothing wrong to be put outside)
                    DynamicList<label> faceVerts(4);

                    label pointi = f[fp];

                    label nextpointi = f[f.fcIndex(fp)];
                    label edgeI = meshTools::findEdge (mesh_, pointi, nextpointi);

                    if (edgeMidPoint[edgeI] >=0)
                    {
                        faceVerts.append(pointi);
                        //change
                        faceVerts.append(edgeMidPoint[edgeI]);

                        //use to check
                        // const edge& e = mesh_.edges()[edgeI];
                        // Pout<< "add:" <<edgeMidPoint[edgeI]<< "is the mid of "<< e <<endl;

                        DynamicList<label> storage;
                        const labelList& fEdges = mesh_.faceEdges(facei,storage);
                        forAll(fEdges,i)
                        {
                            label edgeJ = fEdges[i];
                            if(edgeMidPoint[edgeJ] >= 0 && edgeJ!=edgeI )
                            {
                                faceVerts.append(edgeMidPoint[edgeJ]);

                                //to check
                                // const edge& e = mesh_.edges()[edgeJ];
                                // Pout<< "add:" <<edgeMidPoint[edgeJ]<< "is the mid of "<< e <<endl;
                                break;
                            }
                        }

                        label pointJ = f[f.rcIndex(fp)];
                        
                        faceVerts.append(pointJ);

                        newFace.transfer(faceVerts);

                        //use to check
                        // Pout<< "Split face:" << facei << " verts:" << f
                        //     << " into quad:" << newFace << endl;

                        // Get new owner/neighbour
                        label nei = -1;
                        label own = myGetAnchorCell
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            edgeMidPoint,
                            mesh_.faceOwner()[facei],
                            facei,
                            pointi
                        );

                        if (mesh_.isInternalFace(facei))
                        {
                            nei = myGetAnchorCell
                            (
                                cellAnchorPoints,
                                cellAddedCells,
                                edgeMidPoint,
                                mesh_.faceNeighbour()[facei],
                                facei,
                                pointi
                            );
                        }
                        else
                        {
                            nei = -1;
                        }

                        // use to check
                        // Pout<< "pointi is :" << pointi << endl;
                        // Pout<< "owner cell:" << own << endl;
                        // Pout<< "neighbour cell:" << nei << endl;

                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];

                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];

                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }


                        if (!modifiedFace)
                        {
                            modifiedFace = true;

                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }

            }
            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }

    // use to check
    // Pout<< "done so far"<< endl;

    // 2. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }

    forAll(affectedFace, facei)
    {
        if (affectedFace.get(facei))
        {
            const face& f = mesh_.faces()[facei];
            //use to check
            // Pout<< "affectedFace face:" << facei << " verts:" << f
            // << endl;

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            label anchorFp = findMinLevel(f);
            
            label cello = mesh_.faceOwner()[facei];
            const cell& ocFaces = mesh_.cells()[cello];
            
            label faceowner = -1;
            bool findRightOwnerFace = false;
            forAll(ocFaces, o)
            {
                label faceo = ocFaces[o];

                if(faceMidPoint[faceo] >= 0)
                {
                    const face& of = mesh_.faces()[faceo];
                    forAll(of,fpo)
                    {
                        if( of[fpo] == f[anchorFp] )
                        {
                            findRightOwnerFace = true;
                            faceowner = faceo;
                            break;
                        }
                    }
                }
                if(findRightOwnerFace)
                {
                    break;
                }
            }
            
            //use to check
            // const face& fo = mesh_.faces()[faceowner];
            // Pout<< "chosen faceOwner is :" << fo
            // << endl;

            label nei = -1;
            label own = myGetAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                edgeMidPoint,
                mesh_.faceOwner()[facei],
                faceowner,
                f[anchorFp]
            );

            if (mesh_.isInternalFace(facei))
            {
                label celln = mesh_.faceNeighbour()[facei];
                const cell& ncFaces = mesh_.cells()[celln];
                
                label faceneighbour = -1;
                bool findRightNeighbourFace = false;
                forAll(ncFaces, n)
                {
                    label facen = ncFaces[n];

                    if(faceMidPoint[facen] >= 0)
                    {
                        const face& nf = mesh_.faces()[facen];
                        forAll(nf,fpn)
                        {
                            if( nf[fpn] == f[anchorFp] )
                            {
                                findRightNeighbourFace = true;
                                faceneighbour = facen;
                                break;
                            }
                        }
                    }
                    if(findRightNeighbourFace)
                    {
                        break;
                    }
                }
                
                //use to check
                // const face& fk = mesh_.faces()[faceneighbour];
                // Pout<< "chosen faceNeighbour is :" << fk
                // << endl;

                nei = myGetAnchorCell
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    edgeMidPoint,
                    mesh_.faceNeighbour()[facei],
                    faceneighbour,
                    f[anchorFp]
                );
            }
            else
            {
                nei = -1;
            }

            modFace(meshMod, facei, f, own, nei);

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }


    // 3. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef2::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, celli)
    {
        if (cellMidPoint[celli] >= 0)
        {
            label addPoint = 0;
            label anchorPointi = -1;
            label nextAnchorPointi = -1;
            label lastEdge = -1;
            label lastFace = -1;
            label faceI = -1;
            pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));
            DynamicList<label> newFaceVerts(4);

            bool FirstEdge = true;
            label i = 0;
            const cell& cFaces = mesh_.cells()[celli];
            while(true)
            {
                label facei = cFaces[i];
                //the first divisible edge : add first mid point and find anchor
                if(faceMidPoint[facei] >= 0 && FirstEdge )
                {
                    const face& f = mesh_.faces()[facei];
                    faceI = facei;
                    lastFace = faceI;

                    //decide two anchor point
                    anchorPointi = f[0];
                    nextAnchorPointi = f[f.fcIndex(0)];//start with forward
                    label edgeI = meshTools::findEdge (mesh_, anchorPointi, nextAnchorPointi);

                    if(edgeMidPoint[edgeI] >= 0)
                    {
                        edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
                        newFaceVerts.append(edgeMidPoint[edgeI]);
                        addPoint++;
                        lastEdge = edgeI;
                        FirstEdge = false;
                    }
                    else
                    {
                        nextAnchorPointi = f[f.rcIndex(0)]; //change to backward
                        label edgeI = meshTools::findEdge (mesh_, anchorPointi, nextAnchorPointi);

                        edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
                        newFaceVerts.append(edgeMidPoint[edgeI]);

                        //use to check:
                        // const edge& e = mesh_.edges()[edgeI];
                        // Pout<< "add internal face use edge :" <<e<< endl;
                       
                        if(edgeMidPoint[edgeI] <0 )
                        {
                            Pout<< "something wrong !" << endl;
                        }
                        addPoint++;
                        lastEdge = edgeI;
                        FirstEdge = false;
                    } 
                }
                
                //find the other points to make the new internal face
                if(faceMidPoint[facei] >= 0 && !FirstEdge && facei!=lastFace)
                {
                    const edge& e = mesh_.edges()[lastEdge];
                    const face& f = mesh_.faces()[facei];
                    bool rightFace = false;

                    forAll(f,fp0)
                    {
                        if( f[fp0] == e[0])
                        {
                            if( f[f.fcIndex(fp0)] == e[1] || f[f.rcIndex(fp0)] == e[1])
                            {
                                rightFace = true;
                                break;
                            }
                        }
                    }
                    if(rightFace)
                    {
                        label point0 = -1;
                        label point1 = -1;
                        forAll(f,fp0)
                        {
                            if( f[fp0] != e[0] && f[fp0] != e[1] )
                            {
                                if(point0 == -1)
                                {
                                    point0 = f[fp0];
                                }
                                else
                                {
                                    point1 = f[fp0];
                                    break;
                                }
                            }
                        }
                        label edgeII = meshTools::findEdge (mesh_, point0, point1);

                        edgeMids[edgeII] = mesh_.edges()[edgeII].centre(mesh_.points());
                        newFaceVerts.append(edgeMidPoint[edgeII]);

                        //use to check:
                        // const edge& e = mesh_.edges()[edgeII];
                        // Pout<< "add internal face use edge :" <<e<< endl;
                        
                        if(edgeMidPoint[edgeII] <0 )
                        {
                            Pout<< "something wrong !" << endl;
                        }
                        addPoint++;
                        lastEdge = edgeII;
                        lastFace = facei;
                        FirstEdge = false;
                    }
                }
                i++;
                if(i > 5)
                {
                    i = 0;
                }
                if( addPoint == 4)
                {
                    break;
                }
            }  
            face newFace;
            newFace.transfer(newFaceVerts);
            //use to check
            // Pout<< "done so far 2!,newInternalFace is "<< newFace <<endl;

            label anchorCell0 = myGetAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                edgeMidPoint,
                celli,
                faceI,
                anchorPointi
            );
            //use to check
            // Pout<< "anchorCell0 =" << anchorCell0 << endl;
            label anchorCell1 = myGetAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                edgeMidPoint,
                celli,
                faceI,
                nextAnchorPointi
            );

            //use to check
            // Pout<< "anchorCell1 =" << anchorCell1 << endl;
            // Pout<< "done so far 2.2!" << endl;

            label own, nei;
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                own = anchorCell0;
                nei = anchorCell1;

                ownPt = mesh_.points()[anchorPointi];
                neiPt = mesh_.points()[nextAnchorPointi];

            }
            else
            {
                own = anchorCell1;
                nei = anchorCell0;
                newFace.flip();

                ownPt = mesh_.points()[nextAnchorPointi];
                neiPt = mesh_.points()[anchorPointi];
            }

            if (debug)
            {
                point ownPt, neiPt;

                if (anchorCell0 < anchorCell1)
                {
                    ownPt = mesh_.points()[anchorPointi];
                    neiPt = mesh_.points()[nextAnchorPointi];
                }
                else
                {
                    ownPt = mesh_.points()[nextAnchorPointi];
                    neiPt = mesh_.points()[anchorPointi];
                }

                checkInternalOrientation
                (
                    meshMod,
                    celli,
                    faceI,
                    ownPt,
                    neiPt,
                    newFace
                );
            }

            //use to check
            // Pout<< "done so far 3!"<< endl;

            addInternalFace
            (
                meshMod,
                faceI,
                anchorPointi,
                newFace,
                own,
                nei
            );

        }

    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointi = labelMax;
        label maxPointi = labelMin;

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                minPointi = min(minPointi, cellMidPoint[celli]);
                maxPointi = max(maxPointi, cellMidPoint[celli]);
            }
        }
        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                minPointi = min(minPointi, faceMidPoint[facei]);
                maxPointi = max(maxPointi, faceMidPoint[facei]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointi = min(minPointi, edgeMidPoint[edgeI]);
                maxPointi = max(maxPointi, edgeMidPoint[edgeI]);
            }
        }

        if (minPointi != labelMax && minPointi != mesh_.nPoints())
        {
            FatalErrorInFunction
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointi:" << minPointi
                << " maxPointi:" << maxPointi
                << abort(FatalError);
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    //use to check
    // Pout<< "done setInstance"<< endl;


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexRef2::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, celli)
        {
            const labelList& addedCells = cellAddedCells[celli];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(celli, addedCells);
            }
        }
    }

    //use to check
    // Pout<< "done New unrefinement structure"<< endl;

    // Compact cellAddedCells.

    labelListList refinedCells(cellLabels.size());

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];
        
        refinedCells[i].transfer(cellAddedCells[celli]);
    }

    //use to check
    // Pout<< "done Compact cellAddedCells"<< endl;
    return refinedCells;
}

Foam::labelList Foam::hexRef2::selectUnrefineElems
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitFaces(getSplitElems());

    DynamicList<label> newSplitFaces(splitFaces.size());

    forAll(splitFaces, i)
    {
        label facej = splitFaces[i];

        const face& f = mesh_.faces()[facej];

        forAll(f, i)
        {
            label pointi = f[i];

            bool hasMarked = true;

            if (pFld[pointi] < unrefineLevel)
            {
                // Check that all cells are not marked
                const labelList& pCells = mesh_.pointCells()[pointi];

                hasMarked = false;

                forAll(pCells, pCelli)
                {
                    if (markedCell.get(pCells[pCelli]))
                    {
                        hasMarked = true;
                        break;
                    }
                }
            }

            if (!hasMarked)
            {
                newSplitFaces.append(facej);
                break;
            }
        }
    }


    newSplitFaces.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        consistentUnrefinement
        (
            newSplitFaces,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split faces out of a possible "
        << returnReduce(splitFaces.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}

Foam::labelList Foam::hexRef2::consistentUnrefinement
(
    const labelList& elemsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef2::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorInFunction
            << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // For hexRef2, unrefinement is based on faces
    const labelList& facesToUnrefine(elemsToUnrefine);

    // Loop, modifying facesToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect faces to refine
    // maxSet = true: select faces to refine

    // Maintain boolList for facesToUnrefine and cellsToUnrefine
    PackedBoolList unrefineFace(mesh_.nFaces());

    forAll(facesToUnrefine, i)
    {
        label facei = facesToUnrefine[i];

        unrefineFace.set(facei);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // const labelListList& edgeCells = mesh_.edgeCells();

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefineFace, facei)
        {
            if (unrefineFace.get(facei))
            {
                label pCells1 = mesh_.faceOwner()[facei];
                label pCells2 = mesh_.faceNeighbour()[facei];
                unrefineCell.set(pCells1);
                unrefineCell.set(pCells2);
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            if (ownLevel < (neiLevel-1))
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    // could also combine with unset:
                    // if (!unrefineCell.unset(own))
                    // {
                    //     FatalErrorInFunction
                    //         << "problem cell already unset"
                    //         << abort(FatalError);
                    // }
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if (neiLevel < (ownLevel-1))
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            if (ownLevel < (neiLevel[i]-1))
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if (neiLevel[i] < (ownLevel-1))
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef2::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any edge whose cell neighbour cannot be unrefined.
        forAll(unrefineFace, facei)
        {
            if (unrefineFace.get(facei))
            {
                label pCells1 = mesh_.faceOwner()[facei];
                label pCells2 = mesh_.faceNeighbour()[facei];

                if (!unrefineCell.get(pCells1) || !unrefineCell.get(pCells2) )
                {
                    unrefineFace.unset(facei);
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefineFace, facei)
    {
        if (unrefineFace.get(facei))
        {
            nSet++;
        }
    }

    labelList newFacesToUnrefine(nSet);
    nSet = 0;

    forAll(unrefineFace, facei)
    {
        if (unrefineFace.get(facei))
        {
            newFacesToUnrefine[nSet++] = facei;
        }
    }

    return newFacesToUnrefine;
}

void Foam::hexRef2::calcFaceToSplitPoint
(
    const labelList& splitElems,
    Map<label>& faceToSplitPoint
)
{
    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // For hexRef2, the split elems are faces
    // const labelList& splitPoints(splitElems);
    const labelList& splitFaces(splitElems);

    faceToSplitPoint.resize(splitFaces.size());

    {
        forAll(splitFaces, i)
        {
            label facei = splitFaces[i];

            const labelList& fEdges = mesh_.faceEdges(facei);

            forAll(fEdges, j)
            {
                const edge& e = mesh_.edges()[fEdges[j]];

                forAll(e, k)
                {
                    label pointi = e[k];

                    const labelList& pFaces = mesh_.pointFaces()[pointi];

                    forAll(pFaces, pFacei)
                    {
                        faceToSplitPoint.insert(pFaces[pFacei], pointi);
                    }
                }                
            }
        }
    }
}

// For hexRef2, split elements is face
Foam::labelList Foam::hexRef2::getSplitElems() const
{
    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexRef2::getSplitElems :"
            << " Calculating unrefineable points" << endl;
    }


    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split point
    // >= label of master cell
    labelList splitMaster(mesh_.nFaces(), -1);
    labelList splitMasterLevel(mesh_.nFaces(), 0);

    // Unmark all with not 2 cells
    //const labelListList& pointCells = mesh_.pointCells();
    // maybe useless in hexRef2
    // const labelListList& edgeCells = mesh_.edgeCells();

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, celli)
    {
        const cell& cFaces = mesh_.cells()[celli];

        if (visibleCells[celli] != -1 && history_.parentIndex(celli) >= 0)
        {
            label parentIndex = history_.parentIndex(celli);

            // Check same master.
            forAll(cFaces, i)
            {
                label facei = cFaces[i];

                label masterCelli = splitMaster[facei];

                if (masterCelli == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to celli). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[facei] = parentIndex;
                    splitMasterLevel[facei] = cellLevel_[celli] - 1;
                }
                else if (masterCelli == -2)
                {
                    // Already decided that edge is not splitEdge
                }
                else if
                (
                    (masterCelli != parentIndex)
                 || (splitMasterLevel[facei] != cellLevel_[celli] - 1)
                )
                {
                    // Different masters so edge is on two refinement
                    // patterns
                    splitMaster[facei] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cFaces, i)
            {
                label facei = cFaces[i];

                splitMaster[facei] = -2;
            }
        }
    }

    // Unmark boundary faces
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        splitMaster[facei] = -2;
    }


    // Collect into labelList

    label nSplitFaces = 0;

    forAll(splitMaster, facei)
    {
        if (splitMaster[facei] >= 0)
        {
            nSplitFaces++;
        }
    }

    labelList splitFaces(nSplitFaces);
    nSplitFaces = 0;

    forAll(splitMaster, facei)
    {
        if (splitMaster[facei] >= 0)
        {
            splitFaces[nSplitFaces++] = facei;
        }
    }

    return splitFaces;
}


void Foam::hexRef2::setUnrefinement
(
    const labelList& splitElemLabels,
    polyTopoChange& meshMod
)
{
    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // For hexRef2, unrefinement is based on faces
    const labelList& splitFaceLabels(splitElemLabels);


    //YO- This debug info would require writing edgeSets, which we do not have
    //if (debug)
    //{
    //    Pout<< "hexRef::setUnrefinement :"
    //        << " Checking initial mesh just to make sure" << endl;

    //    checkMesh();

    //    forAll(cellLevel_, celli)
    //    {
    //        if (cellLevel_[celli] < 0)
    //        {
    //            FatalErrorInFunction
    //                << "Illegal cell level " << cellLevel_[celli]
    //                << " for cell " << celli
    //                << abort(FatalError);
    //        }
    //    }


    //    // Write to sets.
    //    pointSet pSet(mesh_, "splitPoints", splitPointLabels);
    //    pSet.write();

    //    cellSet cSet(mesh_, "splitPointCells", splitPointLabels.size());

    //    forAll(splitPointLabels, i)
    //    {
    //        const labelList& pCells = mesh_.pointCells(splitPointLabels[i]);

    //        forAll(pCells, j)
    //        {
    //            cSet.insert(pCells[j]);
    //        }
    //    }
    //    cSet.write();

    //    Pout<< "hexRef::setRefinement : Dumping " << pSet.size()
    //        << " points and "
    //        << cSet.size() << " cells for unrefinement to" << nl
    //        << "    pointSet " << pSet.objectPath() << nl
    //        << "    cellSet " << cSet.objectPath()
    //        << endl;
    //}


    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        labelHashSet splitFaces(splitFaceLabels.size());

        forAll(splitFaceLabels, i)
        {
            splitFaces.insert(splitFaceLabels[i]);
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorInFunction
                << "Initial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitFaceLabels, i)
    {
        label facei = splitFaceLabels[i];

        // Get original cell label

        // const labelList& eCells = mesh_.edgeCells()[edgei];

        // Check
        if ( mesh_.faceNeighbour()[facei] == -1)
        {
            FatalErrorInFunction
                << "splitFace " << facei
                << " should have " << 2
                << " cells using it. It has " << 1
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            // label masterCelli = min(eCells);
            label masterCelli = mesh_.faceOwner()[facei];
            label neighbourCelli = mesh_.faceNeighbour()[facei];

            label regionOwner = cellRegion[masterCelli];
            if (regionOwner == -1)
            {
                FatalErrorInFunction
                    << "Ininitial set of split faces to unrefine does not"
                    << " seem to be consistent or not mid faces"
                    << " of refined cells" << nl
                    << "cell:" << masterCelli << " on splitFace " << facei
                    << " has no region to be merged into"
                    << abort(FatalError);
            }

            if (masterCelli != cellRegionMaster[regionOwner])
            {
                FatalErrorInFunction
                    << "cell:" << masterCelli << " on splitFace:" << facei
                    << " in region " << regionOwner
                    << " has master:" << cellRegionMaster[regionOwner]
                    << " which is not the lowest numbered cell"
                    << abort(FatalError);
            }

            label regionNeighbour = cellRegion[neighbourCelli];
            if (regionNeighbour == -1)
            {
                FatalErrorInFunction
                    << "Ininitial set of split faces to unrefine does not"
                    << " seem to be consistent or not mid faces"
                    << " of refined cells" << nl
                    << "cell:" << neighbourCelli << " on splitFace " << facei
                    << " has no region to be merged into"
                    << abort(FatalError);
            }

            if (masterCelli != cellRegionMaster[regionNeighbour])
            {
                FatalErrorInFunction
                    << "cell:" << neighbourCelli << " on splitFace:" << facei
                    << " in region " << regionNeighbour
                    << " has master:" << cellRegionMaster[regionNeighbour]
                    << " which is not the lowest numbered cell"
                    << abort(FatalError);
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Remove the n cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitFaceLabels, i)
    {
        label facei = splitFaceLabels[i];

        labelList eCells(2);

        label masterCelli = mesh_.faceOwner()[facei];
        eCells[0] = masterCelli;
        label neighbourCelli = mesh_.faceNeighbour()[facei];
        eCells[1] = neighbourCelli;

        cellLevel_[masterCelli]--;
        cellLevel_[neighbourCelli]--;

        history_.combineCells(masterCelli, eCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// ************************************************************************* //
