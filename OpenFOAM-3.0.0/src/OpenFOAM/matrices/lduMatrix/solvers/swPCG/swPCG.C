/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "swPCG.H"
#include "sw_struct.h"
#include "vectorOps.h"
#include "mpi.h"

#include <stdio.h>

#if(SWTIMER)
#include "Timers.hpp"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swPCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<swPCG>
        addswPCGSymMatrixConstructorToTable_;

    amul_translate_array ** swPCG::_matrix_traslate;
    amul_para *swPCG::_amul_parameter;
    int swPCG::if_first = 1;
    refilltion *swPCG::_refill;
    int swPCG::_coarseLevel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swPCG::swPCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{
    register  label nFaces = matrix.upper().size();
    register label nCells = matrix.diag().size();
    const scalar* upperPtr = matrix.upper().begin();
    const scalar* lowerPtr = matrix.lower().begin();
    const label* uPtr = matrix.lduAddr().upperAddr().begin();
    const label* lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* diagPtr =matrix.diag().begin();
    const int coarseLevel_num = 0 ;
    if(if_first)
    {
        _matrix_traslate = (amul_translate_array**) malloc(sizeof(amul_translate_array*)*(coarseLevel_num+1));

        _amul_parameter = (amul_para*)malloc(sizeof(amul_para)*(coarseLevel_num+1));
        _refill = (refilltion*)malloc(sizeof(refilltion)*(coarseLevel_num+1));

        _amul_parameter[0].diagPtr = (scalar*)(&diagPtr[0]);

        _amul_parameter[0].lowerPtr = (scalar*)(&lowerPtr[0]);
        _amul_parameter[0].upperPtr = (scalar*)(&upperPtr[0]);
        _amul_parameter[0].lPtr = (LABEL*)(&lPtr[0]);
        _amul_parameter[0].uPtr = (LABEL*)(&uPtr[0]);
        _amul_parameter[0].nFaces = nFaces;
        _amul_parameter[0].nCells = nCells;

        int _amul_parameter_clumn_size = ROUNDING_UP(nCells, 8);

        _refill[0] = translate_matrix_matrix(_matrix_traslate[0],
                                             &_amul_parameter[0],
                                             _amul_parameter_clumn_size,
                                             8);
        if_first = 0;
    }
    else
    {
        diagPtr = matrix.diag().begin();

        _amul_parameter[0].diagPtr = (scalar*)(&diagPtr[0]);

        _refill[0].lowerPtr = (scalar*)&(lowerPtr[0]);
        _refill[0].upperPtr = (scalar*)&(upperPtr[0]);
        translate_refill_data(&_refill[0]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::swPCG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
#if(SWTIMER)
Timer::startTimer("swPCG");
#endif

    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    label nCells = psi.size();

    scalar* __restrict__ psiPtr = psi.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    scalar wArA = solverPerf.great_;
    scalar wArAold = wArA;

    const scalar* __restrict__ bPtr = source.begin();

    MVM_Arrays arrays1;

    // --- Calculate A.psi
    // matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                            wA,
                            psi,
                            interfaceBouCoeffs_,
                            interfaces_,
                            cmpt,
                            &(_amul_parameter[0]),
                            _matrix_traslate[0],
                            0);

    // --- Calculate initial residual field
    // scalarField rA(source - wA);
    // scalar* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    // scalar normFactor = this->normFactor(psi, source, wA, pA);

    // new
    scalarField rA(nCells);
    scalar* __restrict__ rAPtr = rA.begin();
    scalar vec_temp1[2], vec_temp2[2];
    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            rAPtr[cell]  = bPtr[cell] - wAPtr[cell];
        }
        vec_temp1[0] = sumMag(rA);
        vec_temp1[1] = sumMag(source);
    }
    else
    {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (scalar*)rAPtr;
        arrays1.A2Ptr = (scalar*)bPtr;
        arrays1.A3Ptr = (scalar*)wAPtr;
        arrays1.k1Ptr = &vec_temp1[0];
        arrays1.k2Ptr = &vec_temp1[1];
        // rAPtr = bPtr - wAPtr
        // normFactor += fabs(source)
        residualNormFactor_host(&arrays1);
    }

    if(UPstream::parRun())
    {
        MPI_Allreduce(&vec_temp1, &vec_temp2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else
    {
        vec_temp2[0] = vec_temp1[0];
        vec_temp2[1] = vec_temp1[1];
    }

    const scalar normFactor = 2 * vec_temp2[1] + solverPerf.small_;
    solverPerf.initialResidual() = vec_temp2[0] / normFactor;






    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    // solverPerf.initialResidual() =
    //     gSumMag(rA, matrix().mesh().comm())
    //    /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        // --- Select and construct the preconditioner
        autoPtr<lduMatrix::preconditioner> preconPtr =
        lduMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
#if(SWTIMER)
Timer::startTimer("swPCG preconditioner");
#endif
            preconPtr->precondition(wA, rA, cmpt);

#if(SWTIMER)
Timer::endTimer("swPCG preconditioner");
#endif

            // --- Update search directions:
            // wArA = gSumProd(wA, rA, matrix().mesh().comm());

            wArA = 0.0;
            if(nCells < accUsingSize)
            {
                wArA = gSumProd(wA, rA);
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A2Ptr = (scalar*)wAPtr;
                arrays1.A3Ptr = (scalar*)rAPtr;
                arrays1.k1Ptr = &wArA;
                // wArA += wAPtr * rAPtr
                gSum_host(&arrays1, &slave_userFunc_sumProd);

                reduce(wArA, sumOp<scalar>());
            }

            if (solverPerf.nIterations() == 0)
            {
                // for (label cell=0; cell<nCells; cell++)
                // {
                //     pAPtr[cell] = wAPtr[cell];
                // }

                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell]  = wAPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (scalar*)pAPtr;
                    arrays1.A2Ptr = (scalar*)wAPtr;
                    // pAPtr = wAPtr
                    vectorCopy_host(&arrays1);
                }
            }
            else
            {
                scalar beta = wArA/wArAold;

                // for (label cell=0; cell<nCells; cell++)
                // {
                //     pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                // }


                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (scalar*)pAPtr;
                    arrays1.A2Ptr = (scalar*)wAPtr;
                    arrays1.k1    = beta;
                    // pAPtr = wAPtr + beta*pAPtr
                    vectorOps_host(&arrays1, &slave_userFunc_aEbPk1Mua);
                }
            }


            // --- Update preconditioned residual
            // matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                               wA,
                                               pA,
                                               interfaceBouCoeffs_,
                                               interfaces_,
                                               cmpt,
                                               &(_amul_parameter[0]),
                                               _matrix_traslate[0],
                                               0);


            // scalar wApA = gSumProd(wA, pA, matrix().mesh().comm());
            scalar wApA = 0.0;
            if(nCells < accUsingSize)
            {
                wApA = gSumProd(wA, pA);
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A2Ptr = (scalar*)wAPtr;
                arrays1.A3Ptr = (scalar*)pAPtr;
                arrays1.k1Ptr = &wApA;
                // wApA += wAPtr * pAPtr
                gSum_host(&arrays1, &slave_userFunc_sumProd);

                reduce(wApA, sumOp<scalar>());
            }


            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


            // --- Update solution and residual:

            scalar alpha = wArA/wApA;

            // for (label cell=0; cell<nCells; cell++)
            // {
            //     psiPtr[cell] += alpha*pAPtr[cell];
            //     rAPtr[cell] -= alpha*wAPtr[cell];
            // }

            // solverPerf.finalResidual() =
            //     gSumMag(rA, matrix().mesh().comm())
            //    /normFactor;


            if(nCells < accUsingSize)
            {
                for (label cell=0; cell<nCells; cell++)
                {
                    psiPtr[cell] += alpha*pAPtr[cell];
                    rAPtr[cell]  -= alpha*wAPtr[cell];
                }
                solverPerf.finalResidual() = gSumMag(rA) / normFactor;
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (scalar*)psiPtr;
                arrays1.A2Ptr = (scalar*)pAPtr;
                arrays1.k1    = alpha;
                // psiPtr += k1 * pAPtr
                vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Mub);

                scalar rATemp;
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (scalar*)rAPtr;
                arrays1.A2Ptr = (scalar*)rAPtr;
                arrays1.A3Ptr = (scalar*)wAPtr;
                arrays1.k1Ptr = &rATemp;
                arrays1.k1    = alpha;
                // rAPtr = rAPtr - k1 * wAPtr
                gSum_host(&arrays1, &slave_userFunc_residualSumK);

                reduce(rATemp, sumOp<scalar>());
                solverPerf.finalResidual() = rATemp / normFactor;
            }

        } while
        (
            (
                solverPerf.nIterations()++ < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }
#if(SWTIMER)
Timer::endTimer("swPCG");
#endif
    return solverPerf;
}


// ************************************************************************* //
