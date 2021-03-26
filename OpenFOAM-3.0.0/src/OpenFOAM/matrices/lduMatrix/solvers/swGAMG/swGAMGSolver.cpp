/*---------------------------------------------------------------------------*\
Copyright (C) 2011 OpenFOAM Foundation
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

#include "swGAMGSolver.hpp"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "swAmulTranslate.hpp"
#include "swChebyIfpack2srPCG.hpp"
#include <sstream>
#include <assert.h>
#include "amulMacros.h"
#include "vectorOps.h"
#include "mpi.h"

#include "Timer.hpp"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SW_GAMGSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<SW_GAMGSolver>
        addSW_GAMGSolverMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<SW_GAMGSolver>
        addSW_GAMGAsymSolverMatrixConstructorToTable_;

    amul_translate_array ** SW_GAMGSolver::_matrix_traslate;
    amul_para *SW_GAMGSolver::_amul_parameter;
    int SW_GAMGSolver::if_first = 1;
    refilltion *SW_GAMGSolver::_refill;
    int SW_GAMGSolver::_coarseLevel;

    int* SW_GAMGSolver::upperAddr_int32;
    int* SW_GAMGSolver::lowerAddr_int32;
    int** SW_GAMGSolver::upperAddrLevels_int32;
    int** SW_GAMGSolver::lowerAddrLevels_int32;    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SW_GAMGSolver::SW_GAMGSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    GAMGSolver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    usingSmootherCoarest_(false)
{

    MAX_SW_USING_CORASE_LEVELS = init_MAX_SW_USING_CORASE_LEVELS();
    // init using SmootherCoarsest

    controlDict_.readIfPresent("usingSmootherCoarest", usingSmootherCoarest_);
    // int using_coaseLevels;

    //my_printf_int(MAX_SW_USING_CORASE_LEVELS);
    //controlDict_.readIfPresent("maxSWUingCoarseLevels",MAX_SW_USING_CORASE_LEVELS);
    if(MAX_SW_USING_CORASE_LEVELS > 0 )
    {
        register  label nFaces = matrix.upper().size();
        const scalar* upperPtr = matrix.upper().begin();
        const scalar* lowerPtr = matrix.lower().begin();
        /*const label* uPtr = matrix.lduAddr().upperAddr().begin();
        const label* lPtr = matrix.lduAddr().lowerAddr().begin();*/

// for OpenFOAM int64,  fine level
        const label* uPtr_of = matrix.lduAddr().upperAddr().begin();
        const label* lPtr_of = matrix.lduAddr().lowerAddr().begin();
        int upper_size = matrix.lduAddr().upperAddr().size();
        int lower_size = matrix.lduAddr().lowerAddr().size();
        upperAddr_int32 = new int[upper_size];
        lowerAddr_int32 = new int[lower_size];
        for (int i = 0; i < upper_size; ++i)
        {
            upperAddr_int32[i] = uPtr_of[i]; 
        }
        for (int i = 0; i < lower_size; ++i)
        {
            lowerAddr_int32[i] = lPtr_of[i];
        }
        int* uPtr = upperAddr_int32;
        int* lPtr = lowerAddr_int32;
// end for OpenFOAM int64

        const scalar* diagPtr =matrix.diag().begin();

        int coarseLevel_num = matrixLevels_.size();
        int coarseLevel_i;
        if(if_first)
        {
            SW_ChebyIfpack2srPCG::SET_MAX_SW_USING_CORASE_LEVELS(MAX_SW_USING_CORASE_LEVELS);
            _matrix_traslate = (amul_translate_array**) malloc(sizeof(amul_translate_array*)*(coarseLevel_num+1));
            _amul_parameter = (amul_para*)malloc(sizeof(amul_para)*(coarseLevel_num+1));
            _refill = (refilltion*)malloc(sizeof(refilltion)*(coarseLevel_num+1));
            register label nCells = matrix.diag().size();

            _amul_parameter[0].diagPtr = (SCALAR*)(&diagPtr[0]);

            _amul_parameter[0].lowerPtr = (SCALAR*)(&lowerPtr[0]);
            _amul_parameter[0].upperPtr = (SCALAR*)(&upperPtr[0]);
            _amul_parameter[0].lPtr = (LABEL*)(&lPtr[0]);
            _amul_parameter[0].uPtr = (LABEL*)(&uPtr[0]);
            _amul_parameter[0].nFaces = nFaces;
            _amul_parameter[0].nCells = nCells;

            int _amul_parameter_clumn_size = ROUNDING_UP(nCells,OFFSET_SIZE);

            _refill[0] = translate_matrix_matrix(_matrix_traslate[0],
                                                 &_amul_parameter[0],
                                                 _amul_parameter_clumn_size,
                                                 OFFSET_SIZE);
// for OpenFOAM int64, coarse levels
            //int** SW_GAMGSolver::upperAddrLevels_int32;
            //int** SW_GAMGSolver::lowerAddrLevels_int32;
            upperAddrLevels_int32 = new int*[MAX_SW_USING_CORASE_LEVELS];
            lowerAddrLevels_int32 = new int*[MAX_SW_USING_CORASE_LEVELS];

            for(coarseLevel_i=0; coarseLevel_i<MAX_SW_USING_CORASE_LEVELS-1; coarseLevel_i++)
            {
                const lduMatrix & tmpmatrix = matrixLevels_[coarseLevel_i];
                nCells = tmpmatrix.diag().size();
                nFaces = tmpmatrix.upper().size();
                upperPtr = tmpmatrix.upper().begin();
                lowerPtr = tmpmatrix.lower().begin();
                //uPtr = tmpmatrix.lduAddr().upperAddr().begin();???
                //lPtr = tmpmatrix.lduAddr().lowerAddr().begin();???

// for OpenFOAM int64, coarse levels
                uPtr_of = tmpmatrix.lduAddr().upperAddr().begin();
                lPtr_of = tmpmatrix.lduAddr().lowerAddr().begin();
                upper_size = tmpmatrix.lduAddr().upperAddr().size();
                lower_size = tmpmatrix.lduAddr().lowerAddr().size();
                upperAddrLevels_int32[coarseLevel_i] = new int[upper_size];
                lowerAddrLevels_int32[coarseLevel_i] = new int[lower_size];
                for (int i = 0; i < upper_size; ++i)
                {
                    upperAddrLevels_int32[coarseLevel_i][i] = uPtr_of[i];
                }
                for (int i = 0; i < lower_size; ++i)
                {
                    lowerAddrLevels_int32[coarseLevel_i][i] = lPtr_of[i];
                }
                uPtr = upperAddrLevels_int32[coarseLevel_i];
                lPtr = lowerAddrLevels_int32[coarseLevel_i];
// end for OpenFOAM int64

                diagPtr = tmpmatrix.diag().begin();

                _amul_parameter[coarseLevel_i+1].diagPtr = (SCALAR*)(&diagPtr[0]);

                _amul_parameter[coarseLevel_i+1].lowerPtr = (SCALAR*)(&lowerPtr[0]);
                _amul_parameter[coarseLevel_i+1].upperPtr = (SCALAR*)(&upperPtr[0]);
                _amul_parameter[coarseLevel_i+1].lPtr = (LABEL*)(&lPtr[0]);
                _amul_parameter[coarseLevel_i+1].uPtr = (LABEL*)(&uPtr[0]);
                _amul_parameter[coarseLevel_i+1].nFaces = nFaces;
                _amul_parameter[coarseLevel_i+1].nCells = nCells;
                _amul_parameter_clumn_size = ROUNDING_UP(nCells, OFFSET_SIZE);
                _refill[coarseLevel_i+1] = translate_matrix_matrix((_matrix_traslate[coarseLevel_i+1]),
                                                                   &(_amul_parameter[coarseLevel_i+1]),
                                                                   _amul_parameter_clumn_size,
                                                                   OFFSET_SIZE);
            }
            if_first = 0 ;
        }
        else
        {
            diagPtr = matrix.diag().begin();
            _amul_parameter[0].diagPtr = (SCALAR*)(&diagPtr[0]);

            _refill[0].lowerPtr = (SCALAR*)&(lowerPtr[0]);
            _refill[0].upperPtr = (SCALAR*)&(upperPtr[0]);
            translate_refill_data(&_refill[0]) ;
            for(coarseLevel_i = 0 ; coarseLevel_i < MAX_SW_USING_CORASE_LEVELS - 1; coarseLevel_i++)
            {
                const lduMatrix & tmpmatrix = matrixLevels_[coarseLevel_i];

                diagPtr = tmpmatrix.diag().begin();
                upperPtr = tmpmatrix.upper().begin();
                lowerPtr = tmpmatrix.lower().begin();

                _amul_parameter[coarseLevel_i+1].diagPtr = (SCALAR*)(&diagPtr[0]);

                _refill[coarseLevel_i+1].lowerPtr = (SCALAR*)&(lowerPtr[0]);
                _refill[coarseLevel_i+1].upperPtr = (SCALAR*)&(upperPtr[0]);
                translate_refill_data(&_refill[coarseLevel_i+1]) ;
            }
        }
    }
}

int Foam::SW_GAMGSolver::init_MAX_SW_USING_CORASE_LEVELS()
{
    int i = 0 ;
    const int coarseLevel_num = matrixLevels_.size();
    for (i=0; i<coarseLevel_num; ++i)
    {
        const lduMatrix & tmpmatrix = matrixLevels_[i];
        int nCells = tmpmatrix.diag().size();
        if(nCells < MAXCELLS )
        {
            break;
        }
    }
    return i + 1 ;
}

Foam::solverPerformance Foam::SW_GAMGSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    label xSize = psi.size();
    MVM_Arrays arrays1;

    // Calculate A.psi used to calculate the initial residual
    scalarField Apsi(xSize);

    // matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                            Apsi,
                            psi,
                            interfaceBouCoeffs_,
                            interfaces_,
                            cmpt,
                            &_amul_parameter[0],
                            _matrix_traslate[0],
                            0);

    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    scalarField finestCorrection(xSize);

    // Calculate normalisation factor and initial residual
    scalarField finestResidual(xSize);

    scalar vec_temp1[2], vec_temp2[2];
    if(xSize < accUsingSize)
    {
        finestResidual = source - Apsi;
        vec_temp1[0] = sumMag(finestResidual);
        vec_temp1[1] = sumMag(source);
    }
    else
    {
        init_MVM_Arrays(&arrays1, xSize);
        arrays1.A1Ptr = (SCALAR*)&(finestResidual[0]);
        arrays1.A2Ptr = (SCALAR*)&(source[0]);
        arrays1.A3Ptr = (SCALAR*)&(Apsi[0]);
        arrays1.k1Ptr = &vec_temp1[0];
        arrays1.k2Ptr = &vec_temp1[1];
        // finestResidual = source - Apsi;
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

    // Calculate normalised residual for convergence test
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        // Create coarse grid correction fields
        PtrList<scalarField> coarseCorrFields;

        // Create coarse grid sources
        PtrList<scalarField> coarseSources;

        // Create the smoothers for all levels
        PtrList<lduMatrix::smoother> smoothers;

        // Initialise the above data structures
        initVcycle(coarseCorrFields, coarseSources, smoothers);

        TIMER_START("GAMG vCycle");
        do
        {
            Vcycle
            (
                smoothers,
                psi,
                source,
                Apsi,
                finestCorrection,
                finestResidual,
                coarseCorrFields,
                coarseSources,
                cmpt
            );

            // Calculate finest level residual field
            // matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                    Apsi,
                                    psi,
                                    interfaceBouCoeffs_,
                                    interfaces_,
                                    cmpt,
                                    &_amul_parameter[0],
                                    _matrix_traslate[0],
                                    0);

            if(xSize < accUsingSize)
            {
                finestResidual = source - Apsi;
                solverPerf.finalResidual() = gSumMag(finestResidual)/normFactor;
            }
            else
            {
                scalar rTemp;
                init_MVM_Arrays(&arrays1, xSize);
                arrays1.A1Ptr = (SCALAR*)&(finestResidual[0]);
                arrays1.A2Ptr = (SCALAR*)&(source[0]);
                arrays1.A3Ptr = (SCALAR*)&(Apsi[0]);
                arrays1.k1Ptr = &rTemp;
                // finestResidual = source - Apsi;
                gSum_host(&arrays1, &slave_userFunc_residualSum);
                reduce(rTemp, sumOp<scalar>());
                solverPerf.finalResidual() = rTemp/normFactor;
            }

            if (debug >= 2)
            {
                solverPerf.print(Info.masterStream(matrix().mesh().comm()));
            }
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );

        TIMER_END("GAMG vCycle");
    }

    return solverPerf;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SW_GAMGSolver::~SW_GAMGSolver()
{

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// ************************************************************************* //
