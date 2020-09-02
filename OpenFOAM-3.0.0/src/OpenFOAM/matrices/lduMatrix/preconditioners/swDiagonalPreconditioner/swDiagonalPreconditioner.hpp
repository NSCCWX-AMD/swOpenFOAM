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

Class
    Foam::SW_diagonalPreconditioner

Description
    Diagonal (Jacobi) preconditioner for both symmetric and
    asymmetric matrices.

    The reciprocal of the diagonal is calculated and stored for reuse
    because on most systems '*' is faster than '/'.

SourceFiles
    SW_diagonalPreconditioner.cpp

References

    [1] Templates for the Solution of Linear Systems: Building Blocks
        for Iterative Methods, R. Barrett, M. Barry, T.F. Chan, J. Demmel,
        J. Donato, J. Dongarra, V. Eijkhout, R. Pozo, C. Romine, and
        Van der Vorst, SIAM, 1994, Philadephia, PA, 2nd edition

    [2] Iterative Methods for Sparse Linear Systems, Y. Saad, SIAM, 2003,
        Philadephia, PA, 2nd edition

\*---------------------------------------------------------------------------*/

#ifndef SW_diagonalPreconditioner_H
#define SW_diagonalPreconditioner_H

#include "lduMatrix.H"

namespace Foam
{

class SW_diagonalPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

    //- The reciprocal diagonal
    scalarField rD_;


    // Private Member Functions

    //- Disallow default bitwise copy construct
    SW_diagonalPreconditioner(const SW_diagonalPreconditioner&);

    //- Disallow default bitwise assignment
    void operator=(const SW_diagonalPreconditioner&);

public:

    //- Runtime type information
    TypeName("swDiagonal");

    //- Construct from matrix components and preconditioner solver controls
    SW_diagonalPreconditioner
    (
        const lduMatrix::solver&,
        const dictionary& solverControlsUnused
    );

    virtual ~SW_diagonalPreconditioner()
    {}


    // Member Functions

    //- Return w the preconditioned form of residual r
    virtual void precondition
    (
        scalarField& w,
        const scalarField& r,
        const direction cmpt=0
    ) const;

    //- Return wT the transpose-matrix preconditioned form of residual rT.
    virtual void preconditionT
    (
        scalarField& wT,
        const scalarField& rT,
        const direction cmpt=0
    ) const
    {
        return precondition(wT, rT, cmpt);
    }
};


} // End namespace Foam

#endif

