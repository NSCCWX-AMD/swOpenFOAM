/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "swPBiCGStab.H"
#include "amulMacros.h"
#include "vectorOps.h"
#include "mpi.h"

#include <stdio.h>

#define accUsingSize_Large 1000000

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SW_PBiCGSTAB, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<SW_PBiCGSTAB>
        addSW_PBiCGSTABSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<SW_PBiCGSTAB>
        addSW_PBiCGSTABAsymMatrixConstructorToTable_;

    amul_translate_array ** SW_PBiCGSTAB::_matrix_traslate;
    amul_para *SW_PBiCGSTAB::_amul_parameter;
    int SW_PBiCGSTAB::if_first = 1;
    refilltion *SW_PBiCGSTAB::_refill;
    int SW_PBiCGSTAB::_coarseLevel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SW_PBiCGSTAB::SW_PBiCGSTAB
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

Foam::solverPerformance Foam::SW_PBiCGSTAB::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    const label nCells = psi.size();

    scalar* __restrict__ psiPtr = psi.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField yA(nCells);
    scalar* __restrict__ yAPtr = yA.begin();

    const scalar* __restrict__ bPtr = source.begin();

    scalarField rA(nCells);
    scalar* __restrict__ rAPtr = rA.begin();

    MVM_Arrays arrays1;

    // --- Calculate A.psi
    // matrix_.Amul(yA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                            yA,
                            psi,
                            interfaceBouCoeffs_,
                            interfaces_,
                            cmpt,
                            &(_amul_parameter[0]),
                            _matrix_traslate[0],
                            0);

    // new
    scalar vec_temp1[2], vec_temp2[2];
    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            rAPtr[cell]  = bPtr[cell] - yAPtr[cell];
        }
        vec_temp1[0] = sumMag(rA);
        vec_temp1[1] = sumMag(source);
    }
    else
    {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (scalar*)rAPtr;
        arrays1.A2Ptr = (scalar*)bPtr;
        arrays1.A3Ptr = (scalar*)yAPtr;
        arrays1.k1Ptr = &vec_temp1[0];
        arrays1.k2Ptr = &vec_temp1[1];
        // rAPtr = bPtr - yAPtr
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
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        scalarField AyA(nCells);
        scalar* __restrict__ AyAPtr = AyA.begin();

        scalarField sA(nCells);
        scalar* __restrict__ sAPtr = sA.begin();

        scalarField zA(nCells);
        scalar* __restrict__ zAPtr = zA.begin();

        scalarField tA(nCells);
        scalar* __restrict__ tAPtr = tA.begin();

        // --- Store initial residual
        scalarField rA0(nCells);
        scalar* __restrict__ rA0Ptr = rA0.begin();

        if(nCells < accUsingSize)
        {
            for (register label cell=0; cell<nCells; cell++)
            {
                rA0Ptr[cell]  = rAPtr[cell];
            }
        }
        else
        {
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr = (scalar*)rA0Ptr;
            arrays1.A2Ptr = (scalar*)rAPtr;
            // rA0Ptr = rAPtr
            vectorCopy_host(&arrays1);
        }

        // --- Initial values not used
        scalar rA0rA = 0;
        scalar alpha = 0;
        scalar omega = 0;

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
            // --- Store previous rA0rA
            const scalar rA0rAold = rA0rA;

            rA0rA = 0.0;
            if(nCells < accUsingSize)
            {
                rA0rA = gSumProd(rA0, rA);
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A2Ptr = (scalar*)rA0Ptr;
                arrays1.A3Ptr = (scalar*)rAPtr;
                arrays1.k1Ptr = &rA0rA;
                // rA0rA += rA0Ptr * rAPtr
                gSum_host(&arrays1, &slave_userFunc_sumProd);

                reduce(rA0rA, sumOp<scalar>());
            }

            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(rA0rA)))
            {
                break;
            }

            // --- Update pA
            if (solverPerf.nIterations() == 0)
            {
                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell]  = rAPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (scalar*)pAPtr;
                    arrays1.A2Ptr = (scalar*)rAPtr;
                    // pAPtr = rAPtr
                    vectorCopy_host(&arrays1);
                }
            }
            else
            {
                // --- Test for singularity
                if (solverPerf.checkSingularity(mag(omega)))
                {
                    break;
                }

                const scalar beta = (rA0rA/rA0rAold)*(alpha/omega);

                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell] = rAPtr[cell] + beta*(pAPtr[cell] - omega*AyAPtr[cell]);
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (scalar*)pAPtr;
                    arrays1.A2Ptr = (scalar*)rAPtr;
                    arrays1.A3Ptr = (scalar*)AyAPtr;
                    arrays1.k1    = beta;
                    arrays1.k2    = omega;
                    // pAPtr = rAPtr + beta*(pAPtr - omega*AyAPtr)
                    vectorOps_host(&arrays1, &slave_userFunc_aEbPk1MuSaMik2MucS);
                }
            }

            // --- Precondition pA
            preconPtr->precondition(yA, pA, cmpt);

            // --- Calculate AyA
            // matrix_.Amul(AyA, yA, interfaceBouCoeffs_, interfaces_, cmpt);
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                    AyA,
                                    yA,
                                    interfaceBouCoeffs_,
                                    interfaces_,
                                    cmpt,
                                    &(_amul_parameter[0]),
                                    _matrix_traslate[0],
                                    0);


            scalar rA0AyA = 0.0;
            if(nCells < accUsingSize)
            {
                rA0AyA = gSumProd(rA0, AyA);
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A2Ptr = (scalar*)rA0Ptr;
                arrays1.A3Ptr = (scalar*)AyAPtr;
                arrays1.k1Ptr = &rA0AyA;
                // rA0AyA += rA0Ptr * AyAPtr
                gSum_host(&arrays1, &slave_userFunc_sumProd);

                reduce(rA0AyA, sumOp<scalar>());
            }

            alpha = rA0rA/rA0AyA;

            // --- Calculate sA
            // --- Test sA for convergence
            if(nCells < accUsingSize)
            {
                for (label cell=0; cell<nCells; cell++)
                {
                    sAPtr[cell] = rAPtr[cell] - alpha*AyAPtr[cell];
                }
                solverPerf.finalResidual() = gSumMag(sA) / normFactor;
            }
            else
            {
                scalar sATemp;
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (scalar*)sAPtr;
                arrays1.A2Ptr = (scalar*)rAPtr;
                arrays1.A3Ptr = (scalar*)AyAPtr;
                arrays1.k1Ptr = &sATemp;
                arrays1.k1    = alpha;
                gSum_host(&arrays1, &slave_userFunc_residualSumK);

                reduce(sATemp, sumOp<scalar>());

                solverPerf.finalResidual() = sATemp / normFactor;
            }

            if (solverPerf.checkConvergence(tolerance_, relTol_))
            {
                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        psiPtr[cell] += alpha*yAPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (scalar*)psiPtr;
                    arrays1.A2Ptr = (scalar*)yAPtr;
                    arrays1.k1    = alpha;
                    // psiPtr += k1 * yAPtr
                    vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Mub);
                }

                solverPerf.nIterations()++;

                return solverPerf;
            }

            // --- Precondition sA
            preconPtr->precondition(zA, sA, cmpt);

            // --- Calculate tA
            // matrix_.Amul(tA, zA, interfaceBouCoeffs_, interfaces_, cmpt);
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                    tA,
                                    zA,
                                    interfaceBouCoeffs_,
                                    interfaces_,
                                    cmpt,
                                    &(_amul_parameter[0]),
                                    _matrix_traslate[0],
                                    0);

            // scalarField tA_tmp(nCells);
            // scalar* __restrict__ tAPtr_tmp = tA_tmp.begin();
            // SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
            //                         tA_tmp,
            //                         zA,
            //                         interfaceBouCoeffs_,
            //                         interfaces_,
            //                         cmpt,
            //                         &(_amul_parameter[0]),
            //                         _matrix_traslate[0],
            //                         0);

            // for(int igu=0; igu<nCells; igu++)
            // {
            //     double diff = fabs((tAPtr[igu] - tAPtr_tmp[igu])/(tAPtr[igu]+solverPerf.small_));
            //     if(diff > 0.000000000001)
            //     {
            //         std::cout << "tA = " << tAPtr[igu] << ", tA_tmp = " <<  tAPtr_tmp[igu] << std::endl;
            //         int myi = 0;
            //         while (0 == myi)
            //         {
            //             Foam::sleep(5);
            //         }
            //     }
            // }

            // for(int igu=0; igu<nCells; igu++)
            // {
            //     tAPtr[igu] = tAPtr_tmp[igu];
            // }



            scalar tAtA = 0.0;
            if(nCells < accUsingSize)
            {
                tAtA = gSumSqr(tA);
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A2Ptr = (scalar*)tAPtr;
                // arrays1.A3Ptr = (scalar*)tAPtr;
                arrays1.k1Ptr = &tAtA;
                // tAtA += tAPtr * tAPtr
                gSum_host(&arrays1, &slave_userFunc_sumSqr);

                reduce(tAtA, sumOp<scalar>());
            }

            // --- Calculate omega from tA and sA
            //     (cheaper than using zA with preconditioned tA)
            if(nCells < accUsingSize)
            {
                omega = gSumProd(tA, sA)/tAtA;
            }
            else
            {
                omega = 0.0;
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A2Ptr = (scalar*)tAPtr;
                arrays1.A3Ptr = (scalar*)sAPtr;
                arrays1.k1Ptr = &omega;
                // omega += tAPtr * sAPtr
                gSum_host(&arrays1, &slave_userFunc_sumProd);

                reduce(omega, sumOp<scalar>());

                omega = omega / tAtA;
            }

            // --- Update solution and residual
            if(nCells < accUsingSize)
            {
                for (label cell=0; cell<nCells; cell++)
                {
                    psiPtr[cell] += alpha*yAPtr[cell] + omega*zAPtr[cell];
                    rAPtr[cell] = sAPtr[cell] - omega*tAPtr[cell];
                }
                solverPerf.finalResidual() = gSumMag(rA) / normFactor;
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (scalar*)psiPtr;
                arrays1.A2Ptr = (scalar*)yAPtr;
                arrays1.A3Ptr = (scalar*)zAPtr;
                arrays1.k1    = alpha;
                arrays1.k2    = omega;
                // psiPtr += k1 * yAPtr + k2 * zAPtr
                vectorOps_host(&arrays1, &slave_userFunc_aEaPk1MubPk2Muc);

                scalar rATemp;
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (scalar*)rAPtr;
                arrays1.A2Ptr = (scalar*)sAPtr;
                arrays1.A3Ptr = (scalar*)tAPtr;
                arrays1.k1Ptr = &rATemp;
                arrays1.k1    = omega;
                // rAPtr = sAPtr - k1 * tAPtr
                gSum_host(&arrays1, &slave_userFunc_residualSumK);

                reduce(rATemp, sumOp<scalar>());

                solverPerf.finalResidual() = rATemp / normFactor;
            }
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    return solverPerf;
}

Foam::SW_PBiCGSTAB::~SW_PBiCGSTAB()
{
}

// ************************************************************************* //
