/*---------------------------------------------------------------------------*\
Copyright (C) 2011 OpenFOAM Foundation
Copyright (C) 2014 Applied CCM
-------------------------------------------------------------------------------
License
    This file is part of CAELUS.

    CAELUS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CAELUS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CAELUS.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "swDILUPreconditioner.hpp"
#include "amulMacros.h"
#include "vectorOps.h"
#include "mpi.h"

namespace Foam
{
    defineTypeNameAndDebug(SW_DILUPreconditioner, 0);

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<SW_DILUPreconditioner>
        addSW_DILUPreconditionerAsymMatrixConstructorToTable_;
}

Foam::SW_DILUPreconditioner::SW_DILUPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag())
{
    calcReciprocalD(rD_, sol.matrix());
}

void Foam::SW_DILUPreconditioner::calcReciprocalD
(
    scalarField& rD,
    const lduMatrix& matrix
)
{
    scalar* __restrict__ rDPtr = rD.begin();

    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = matrix.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix.lower().begin();

    register label nFaces = matrix.upper().size();
    for (register label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= upperPtr[face]*lowerPtr[face]/rDPtr[lPtr[face]];
    }


    // Calculate the reciprocal of the preconditioned diagonal
    register label nCells = rD.size();


    for (register label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }
}

void Foam::SW_DILUPreconditioner::precondition
(
    scalarField& w,
    const scalarField& r,
    const direction
) const
{
    scalar* __restrict__ wPtr = w.begin();
    const scalar* __restrict__ rPtr = r.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().begin();
    const label* const __restrict__ losortPtr =
        solver_.matrix().lduAddr().losortAddr().begin();

    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().begin();
    const scalar* const __restrict__ lowerPtr =
        solver_.matrix().lower().begin();

    register label nCells = w.size();
    register label nFaces = solver_.matrix().upper().size();
    register label nFacesM1 = nFaces - 1;

    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            wPtr[cell] = rDPtr[cell] * rPtr[cell];
        }
    }
    else
    {
        MVM_Arrays arrays1;
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (SCALAR*)wPtr;
        arrays1.A2Ptr = (SCALAR*)rDPtr;
        arrays1.A3Ptr = (SCALAR*)rPtr;
        // wPtr = rDPtr * rPtr
        vectorOps_host(&arrays1, &slave_userFunc_aEbMuc);
    }

    register label sface;

    for (register label face=0; face<nFaces; face++)
    {
        sface = losortPtr[face];
        wPtr[uPtr[sface]] -=
            rDPtr[uPtr[sface]]*lowerPtr[sface]*wPtr[lPtr[sface]];
    }

    for (register label face=nFacesM1; face>=0; face--)
    {
        wPtr[lPtr[face]] -=
            rDPtr[lPtr[face]]*upperPtr[face]*wPtr[uPtr[face]];
    }
}

void Foam::SW_DILUPreconditioner::preconditionT
(
    scalarField& wT,
    const scalarField& rT,
    const direction
) const
{
    scalar* __restrict__ wTPtr = wT.begin();
    const scalar* __restrict__ rTPtr = rT.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().begin();
    const label* const __restrict__ losortPtr =
        solver_.matrix().lduAddr().losortAddr().begin();

    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().begin();
    const scalar* const __restrict__ lowerPtr =
        solver_.matrix().lower().begin();

    register label nCells = wT.size();
    register label nFaces = solver_.matrix().upper().size();
    register label nFacesM1 = nFaces - 1;

    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            wTPtr[cell] = rDPtr[cell]*rTPtr[cell];
        }
    }
    else
    {
        MVM_Arrays arrays1;
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (SCALAR*)wTPtr;
        arrays1.A2Ptr = (SCALAR*)rDPtr;
        arrays1.A3Ptr = (SCALAR*)rTPtr;
        // wTPtr = rDPtr * rTPtr
        vectorOps_host(&arrays1, &slave_userFunc_aEbMuc);
    }

    for (register label face=0; face<nFaces; face++)
    {
        wTPtr[uPtr[face]] -=
            rDPtr[uPtr[face]]*upperPtr[face]*wTPtr[lPtr[face]];
    }

    register label sface;

    for (register label face=nFacesM1; face>=0; face--)
    {
        sface = losortPtr[face];
        wTPtr[lPtr[sface]] -=
            rDPtr[lPtr[sface]]*lowerPtr[sface]*wTPtr[uPtr[sface]];
    }
}

