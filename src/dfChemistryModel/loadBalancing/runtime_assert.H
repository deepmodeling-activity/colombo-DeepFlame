/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland 
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.
    https://github.com/blttkgl/DLBFoam

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

Description
    A runtime assert that does nothing on release builds

\*---------------------------------------------------------------------------*/
#ifndef runtime_assert_H
#define runtime_assert_H

#include <exception>      // std::runtime_error
#include "Ostream.H"

namespace Foam
{

#ifdef DEBUG
static inline void
runtime_assert(bool condition, const char* msg)
{
    if(!condition)
        throw std::runtime_error(msg);
}
#else
static inline void
runtime_assert(
    [[maybe_unused]] bool condition, [[maybe_unused]] const char* msg)
{
}
#endif



} // namespace Foam

#endif