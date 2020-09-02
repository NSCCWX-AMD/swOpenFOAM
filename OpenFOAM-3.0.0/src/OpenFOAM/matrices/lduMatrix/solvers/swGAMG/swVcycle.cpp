#include "swGAMGSolver.hpp"
#include <stdio.h>
#include <iostream>
#include "SubField.H"
#include "amulMacros.h"
#include "vectorOps.h"
#include "BICCG.H"
#include "srPCG.hpp"

void Foam::SW_GAMGSolver::Vcycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    scalarField& psi,
    const scalarField& source,
    scalarField& Apsi,
    scalarField& finestCorrection,
    scalarField& finestResidual,
    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources,
    const direction cmpt
) const
{
    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }

    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        // If the optional pre-smoothing sweeps are selected
        // smooth the coarse-grid field for the restriced source

        if (nPreSweeps_)
        {
            coarseCorrFields[leveli] = 0.0;

            smoothers[leveli + 1].smooth
            (
                coarseCorrFields[leveli],
                coarseSources[leveli],
                cmpt,
                nPreSweeps_ + leveli
            );

            scalarField::subField ACf
            (
                Apsi,
                coarseCorrFields[leveli].size()
            );

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if (scaleCorrection_ && leveli < coarsestLevel - 1)
            {
                scalar sf = scalingFactor
                (
                    const_cast<scalarField&>(ACf.operator const scalarField&()),
                    matrixLevels_[leveli],
                    coarseCorrFields[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources[leveli],
                    cmpt,
                    leveli+1
                );

                if (debug >= 2)
                {
                    Pout<< sf << " ";
                }

                label nCells = coarseCorrFields[leveli].size();
                if(nCells < accUsingSize)
                {
                    coarseCorrFields[leveli] *= sf;
                }
                else
                {
                    MVM_Arrays arrays1;
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (SCALAR*)(coarseCorrFields[leveli].begin());
                    arrays1.k1    = sf;
                    vectorOps_host(&arrays1, &slave_userFunc_aEk1Mua);
                }
            }

            // Correct the residual with the new solution
            // matrixLevels_[leveli].Amul
            // (
            //     const_cast<scalarField&>(ACf.operator const scalarField&()),
            //     coarseCorrFields[leveli],
            //     interfaceLevelsBouCoeffs_[leveli],
            //     interfaceLevels_[leveli],
            //     cmpt
            // );
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrixLevels_[leveli],
                                    const_cast<scalarField&>(ACf.operator const scalarField&()),
                                    coarseCorrFields[leveli],
                                    interfaceLevelsBouCoeffs_[leveli],
                                    interfaceLevels_[leveli],
                                    cmpt,
                                    &_amul_parameter[leveli+1],
                                    _matrix_traslate[leveli+1],
                                    leveli+1);
            coarseSources[leveli] -= ACf;
        }

        // Residual is equal to source
        agglomeration_.restrictField
        (
            coarseSources[leveli + 1],
            coarseSources[leveli],
            leveli + 1,
            true
        );
    }

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< endl;
    }

////////////////// author : hf
// Solve Coarsest level with either an iterative or direct solver

    //label nCellsLevels = coarseCorrFields[coarsestLevel].size();

    // if(usingSmootherCoarest_)
    // {
    //     smoothCoarsestLevel
    //     (
    //         coarseCorrFields[coarsestLevel],
    //         coarseSources[coarsestLevel],
    //         smoothers,
    //         cmpt
    //     );
    // }
    // else
    // {
    //     solveCoarsestLevel
    //     (
    //         coarseCorrFields[coarsestLevel],
    //         coarseSources[coarsestLevel]
    //     );
    // }

    // Solve Coarsest level with either an iterative or direct solver
    SW_solveCoarsestLevel
    (
        coarseCorrFields[coarsestLevel],
        coarseSources[coarsestLevel]
    );

    if (debug >= 2)
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)
    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        label nCells = coarseCorrFields[leveli].size();
        MVM_Arrays arrays1;

        // Create a field for the pre-smoothed correction field
        // as a sub-field of the finestCorrection which is not
        // currently being used
        scalarField::subField preSmoothedCoarseCorrField
        (
            finestCorrection,
            coarseCorrFields[leveli].size()
        );

        // Only store the preSmoothedCoarseCorrField is pre-smoothing is used
        if (nPreSweeps_)
        {
            if(nCells < accUsingSize)
            {
                preSmoothedCoarseCorrField.assign(coarseCorrFields[leveli]);
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (SCALAR*)(preSmoothedCoarseCorrField.begin());
                arrays1.A2Ptr = (SCALAR*)(coarseCorrFields[leveli].begin());
                vectorCopy_host(&arrays1);
            }
        }

        agglomeration_.prolongField
        (
            coarseCorrFields[leveli],
            coarseCorrFields[leveli + 1],
            leveli + 1,
            true
        );

        // Scale coarse-grid correction field
        // but not on the coarsest level because it evaluates to 1
        if (scaleCorrection_ && leveli < coarsestLevel - 1)
        {
            // Create A.psi for this coarse level as a sub-field of Apsi
            scalarField::subField ACf
            (
                Apsi,
                coarseCorrFields[leveli].size()
            );

            scalar sf = scalingFactor
            (
                const_cast<scalarField&>(ACf.operator const scalarField&()),
                matrixLevels_[leveli],
                coarseCorrFields[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                coarseSources[leveli],
                cmpt,
                leveli+1
            );


            if (debug >= 2)
            {
                Pout<< sf << " ";
            }

            if(nCells < accUsingSize)
            {
                coarseCorrFields[leveli] *= sf;
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (SCALAR*)(coarseCorrFields[leveli].begin());
                arrays1.k1    = sf;
                vectorOps_host(&arrays1, &slave_userFunc_aEk1Mua);
            }
        }

        // Only add the preSmoothedCoarseCorrField is pre-smoothing is used
        if (nPreSweeps_)
        {
            if(nCells < accUsingSize)
            {
                coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (SCALAR*)(coarseCorrFields[leveli].begin());
                arrays1.A2Ptr = (SCALAR*)(preSmoothedCoarseCorrField.begin());
                vectorOps_host(&arrays1, &slave_userFunc_aEaPb);
            }
        }

        smoothers[leveli + 1].smooth
        (
            coarseCorrFields[leveli],
            coarseSources[leveli],
            cmpt,
            nPostSweeps_ + leveli
        );
    }

    // Prolong the finest level correction
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true
    );

    if (scaleCorrection_)
    {
        // Calculate finest level scaling factor
        scalar fsf = scalingFactor
        (
            Apsi,
            matrix_,
            finestCorrection,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt,
            0
        );

        if (debug >= 2)
        {
            Pout<< fsf << endl;
        }


        label nCells = psi.size();
        if(nCells < accUsingSize)
        {
            for (register label cell=0; cell<nCells; cell++)
            {
                psi[cell] += fsf*finestCorrection[cell];
            }
        }
        else
        {
            MVM_Arrays arrays1;
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr = (SCALAR*)&(psi[0]);
            arrays1.A2Ptr = (SCALAR*)&(finestCorrection[0]);
            arrays1.k1    = fsf;
            // psi += fsf*finestCorrection;
            vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Mub);
        }
    }
    else
    {
        label nCells = psi.size();
        if(nCells < accUsingSize)
        {
            for (register label cell=0; cell<nCells; cell++)
            {
                psi[cell] += finestCorrection[cell];
            }
        }
        else
        {
            MVM_Arrays arrays1;
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr = (SCALAR*)&(psi[0]);
            arrays1.A2Ptr = (SCALAR*)&(finestCorrection[0]);
            // psi += finestCorrection;
            vectorOps_host(&arrays1, &slave_userFunc_aEaPb);
        }
    }

    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
}

void Foam::SW_GAMGSolver::initVcycle
(
    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources,
    PtrList<lduMatrix::smoother>& smoothers
) const
{
    coarseCorrFields.setSize(matrixLevels_.size());
    coarseSources.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    // Create the smoother for the finest level

    ////writen by Changxi , for the constructor of Cheby to know the coarse level

    _coarseLevel = 0 ;
    smoothers.set
    (
        0,
        lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        )
    );

    forAll(matrixLevels_, leveli)
    {
        coarseCorrFields.set
        (
            leveli,
            new scalarField
            (
                agglomeration_.meshLevel(leveli + 1).lduAddr().size()
            )
        );

        coarseSources.set
        (
            leveli,
            new scalarField
            (
                agglomeration_.meshLevel(leveli + 1).lduAddr().size()
            )
        );
        ////writen by Changxi , for the constructor of Cheby to know the coarse level
        _coarseLevel = leveli + 1;
        ////
        smoothers.set
        (
            leveli + 1,
            lduMatrix::smoother::New
            (
                fieldName_,
                matrixLevels_[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevelsIntCoeffs_[leveli],
                interfaceLevels_[leveli],
                controlDict_
            )
        );
    }
}

void Foam::SW_GAMGSolver::SW_solveCoarsestLevel
(
    scalarField& coarsestCorrField,
    const scalarField& coarsestSource
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;

    label coarseComm = matrixLevels_[coarsestLevel].mesh().comm();
    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = coarseComm;

    if (directSolveCoarsest_)
    {
        coarsestCorrField = coarsestSource;
        coarsestLUMatrixPtr_->solve(coarsestCorrField);
    }
    else
    {
        coarsestCorrField = 0;
        solverPerformance coarseSolverPerf;

        if (matrixLevels_[coarsestLevel].asymmetric())
        {
            coarseSolverPerf = BICCG
            (
                "coarsestLevelCorr",
                matrixLevels_[coarsestLevel],
                interfaceLevelsBouCoeffs_[coarsestLevel],
                interfaceLevelsIntCoeffs_[coarsestLevel],
                interfaceLevels_[coarsestLevel],
                tolerance_,
                relTol_
            ).solve
            (
                coarsestCorrField,
                coarsestSource
            );
        }
        else
        {
            coarseSolverPerf = srPCG
            (
                "coarsestLevelCorr",
                matrixLevels_[coarsestLevel],
                interfaceLevelsBouCoeffs_[coarsestLevel],
                interfaceLevelsIntCoeffs_[coarsestLevel],
                interfaceLevels_[coarsestLevel],
                tolerance_,
                relTol_
            ).solve
            (
                coarsestCorrField,
                coarsestSource
            );
        }

        if (debug >= 2)
        {
            coarseSolverPerf.print(Info.masterStream(coarseComm));
        }
    }

    UPstream::warnComm = oldWarn;
}
