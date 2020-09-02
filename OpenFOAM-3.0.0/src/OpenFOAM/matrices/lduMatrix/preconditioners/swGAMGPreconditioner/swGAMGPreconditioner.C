/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "swGAMGPreconditioner.H"
#include "swAmulTranslate.hpp"
#include "sw_struct.h"
#include "vectorOps.h"
#include "mpi.h"

#if(SWTIMER)
#include "Timers.hpp"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swGAMGPreconditioner, 0);

    lduMatrix::preconditioner::addsymMatrixConstructorToTable
    <swGAMGPreconditioner> addswGAMGPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::addasymMatrixConstructorToTable
    <swGAMGPreconditioner> addswGAMGPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swGAMGPreconditioner::swGAMGPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& solverControls
)
:
    SW_GAMGSolver
    (
        sol.fieldName(),
        sol.matrix(),
        sol.interfaceBouCoeffs(),
        sol.interfaceIntCoeffs(),
        sol.interfaces(),
        solverControls
    ),
    lduMatrix::preconditioner(sol),
    nVcycles_(2)
{
    readControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swGAMGPreconditioner::~swGAMGPreconditioner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swGAMGPreconditioner::readControls()
{
    SW_GAMGSolver::readControls();
    nVcycles_ = controlDict_.lookupOrDefault<label>("nVcycles", 2);
}


void Foam::swGAMGPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction cmpt
) const
{
    wA = 0.0;
    scalarField AwA(wA.size());
    scalarField finestCorrection(wA.size());
    scalarField finestResidual(rA);

    // Create coarse grid correction fields
    PtrList<scalarField> coarseCorrFields;

    // Create coarse grid sources
    PtrList<scalarField> coarseSources;

    // Create the smoothers for all levels
    PtrList<lduMatrix::smoother> smoothers;

    // Scratch fields if processor-agglomerated coarse level meshes
    // are bigger than original. Usually not needed
    // scalarField ApsiScratch;
    // scalarField finestCorrectionScratch;

    // Initialise the above data structures
    initVcycle
    (
        coarseCorrFields,
        coarseSources,
        smoothers
    );

    for (label cycle=0; cycle<nVcycles_; cycle++)
    {
        Vcycle
        (
            smoothers,
            wA,
            rA,
            AwA,
            finestCorrection,
            finestResidual,
            coarseCorrFields,
            coarseSources,
            cmpt
        );

        if (cycle < nVcycles_-1)
        {
            // Calculate finest level residual field
            // matrix_.Amul(AwA, wA, interfaceBouCoeffs_, interfaces_, cmpt);
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                               AwA,
                                               wA,
                                               interfaceBouCoeffs_,
                                               interfaces_,
                                               cmpt,
                                               &_amul_parameter[0],
                                               _matrix_traslate[0],
                                               0);

            const label nCells = rA.size();
            if(nCells < accUsingSize)
            {
                finestResidual = rA - AwA;
            }
            else
            {
                MVM_Arrays arrays1;
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (scalar*)finestResidual.begin();
                arrays1.A2Ptr = (scalar*)rA.begin();
                arrays1.A3Ptr = (scalar*)AwA.begin();
                // finestResidual = rA - AwA
                vectorOps_host(&arrays1, &slave_userFunc_aEbMic);
            }
            // finestResidual = rA;
            // finestResidual -= AwA;
        }
    }
}


// ************************************************************************* //
