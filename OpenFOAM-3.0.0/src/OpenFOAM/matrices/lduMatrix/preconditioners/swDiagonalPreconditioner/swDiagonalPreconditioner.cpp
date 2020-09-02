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

#include "swDiagonalPreconditioner.hpp"
#include "amulMacros.h"
#include "vectorOps.h"

#include "mpi.h"

namespace Foam
{
    defineTypeNameAndDebug(SW_diagonalPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<SW_diagonalPreconditioner>
        addSW_diagonalPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<SW_diagonalPreconditioner>
        addSW_diagonalPreconditionerAsymMatrixConstructorToTable_;
}

Foam::SW_diagonalPreconditioner::SW_diagonalPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag().size())
{
    scalar* __restrict__ rDPtr = rD_.begin();
    const scalar* __restrict__ DPtr = solver_.matrix().diag().begin();

    register label nCells = rD_.size();

    // Generate reciprocal diagonal

    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            rDPtr[cell] = 1.0 / DPtr[cell];
        }
    }
    else
    {
        MVM_Arrays arrays1;
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (SCALAR*)rDPtr;
        arrays1.A2Ptr = (SCALAR*)DPtr;
        // rDPtr = 1 / DPtr
        vectorOps_host(&arrays1, &slave_userFunc_aE1Db);
    }
}

void Foam::SW_diagonalPreconditioner::precondition
(
    scalarField& w,
    const scalarField& r,
    const direction
) const
{
    scalar* __restrict__ wPtr = w.begin();
    const scalar* __restrict__ rPtr = r.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    register label nCells = w.size();

    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            wPtr[cell] = rDPtr[cell]*rPtr[cell];
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
}


