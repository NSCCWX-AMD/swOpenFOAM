#include "swChebyIfpack2srPCG.hpp"
#include "smoothSolver.H"
#include <time.h>
#include <vector>
#include "vector2D.H"
#include "swGAMGSolver.hpp"
#include "amulMacros.h"
#include "vectorOps.h"
#include "mpi.h"

#define accUsingSize_Large 1000000
#define mpiTest 0

namespace Foam
{
    defineTypeNameAndDebug(SW_ChebyIfpack2srPCG, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<SW_ChebyIfpack2srPCG>
        addSW_ChebyIfpack2srPCGSymMatrixConstructorToTable_;
    int SW_ChebyIfpack2srPCG::MAX_SW_USING_CORASE_LEVELS = 0;
}

void  Foam::SW_ChebyIfpack2srPCG::SET_MAX_SW_USING_CORASE_LEVELS(int level)
{
        MAX_SW_USING_CORASE_LEVELS = level;
}

Foam::SW_ChebyIfpack2srPCG::SW_ChebyIfpack2srPCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    ChebyshevSmootherIfpack2
    (
            fieldName,
            matrix,
            interfaceBouCoeffs,
            interfaceIntCoeffs,
            interfaces
    )
{
    _coarseLevel = SW_GAMGSolver::_coarseLevel;

    _amul_parameter = & (SW_GAMGSolver::_amul_parameter[_coarseLevel]);
    _matrix_translate =  (SW_GAMGSolver::_matrix_traslate[_coarseLevel]);

}

Foam::SW_ChebyIfpack2srPCG::~SW_ChebyIfpack2srPCG()
{

}

void Foam::SW_ChebyIfpack2srPCG::smooth
(
    const word& fieldName_,
    scalarField& x,
    const lduMatrix& matrix_,
    const scalarField& b,
    const FieldField<Field, scalar>& interfaceBouCoeffs_,
    const FieldField<Field, scalar>& interfaceIntCoeffs_,
    const lduInterfaceFieldPtrsList& interfaces_,
    const direction cmpt,
    const label nSweeps
) const
{
    if (matrix_.asymmetric())
    {
        Info << "ERROR!!" << endl
             << "Matrix asymmetric!!!" << endl
             << "Cheby smoother cannot work!!!" << endl;
        std::exit(0);
    }

    register const label nCells = x.size();
    scalar* __restrict__ xPtr = x.begin();

    scalarField p(nCells);
    scalar* __restrict__ pPtr = p.begin();

    scalarField w(nCells);
    scalar* __restrict__ wPtr = w.begin();

    scalarField r(nCells);
    scalar* __restrict__ rPtr = r.begin();

    scalarField rD(nCells);
    scalar* __restrict__ rDPtr = rD.begin();

    const scalar* __restrict__ bPtr = b.begin();
    const scalar* __restrict__ DPtr = matrix_.diag().begin();

    scalarField niu(nCells);
    scalar* __restrict__ niuPtr = niu.begin();

    scalarField s(nCells);
    scalar* __restrict__ sPtr = s.begin();

    MVM_Arrays arrays1;

    // diagonal precondition
    if(nCells < accUsingSize)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            rDPtr[cell] = 1.0 / DPtr[cell];
        }
    }
    else
    {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (SCALAR*)rDPtr;
        arrays1.A2Ptr = (SCALAR*)DPtr;
        // rDPtr = 1 / DPtr
        vectorOps_host(&arrays1, &slave_userFunc_aE1Db);
    }

    // compute eigenvalue using diagnoal conjugate gradient sweeps
    if(eigFirstTimeComputed_)
    {
        eigFirstTimeComputed_ = false;

        std::vector<scalar> alphas;
        std::vector<scalar> betas;

        scalar rho = 1e20;
        scalar rhoOld = rho;
        scalar wApA = 0.0;

        // matrix_.Amul(w, x, interfaceBouCoeffs_, interfaces_, cmpt);
        SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                w,
                                x,
                                interfaceBouCoeffs_,
                                interfaces_,
                                cmpt,
                                _amul_parameter,
                                _matrix_translate,
                                _coarseLevel);

        if(nCells < accUsingSize)
        {
            for (register label cell=0; cell<nCells; cell++)
            {
                rPtr[cell] = bPtr[cell] - wPtr[cell];
            }
        }
        else
        {
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr = (SCALAR*)rPtr;
            arrays1.A2Ptr = (SCALAR*)bPtr;
            arrays1.A3Ptr = (SCALAR*)wPtr;
            // rPtr = bPtr - wPtr
            vectorOps_host(&arrays1, &slave_userFunc_aEbMic);
        }

        for(label nIter=0; nIter<nPCGsCheby_; nIter++)
        {
            // Store previous rho
            rhoOld = rho;
            rho    = 0.0;

            // Apply preconditioner
            // Digonal preconditioner and update rho
            if(nCells < accUsingSize)
            {

                for (register label cell=0; cell<nCells; cell++)
                {
                    wPtr[cell] = rDPtr[cell] * rPtr[cell];
                    rho       +=  wPtr[cell] * rPtr[cell];
                }
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (SCALAR*)wPtr;
                arrays1.A2Ptr = (SCALAR*)rDPtr;
                arrays1.A3Ptr = (SCALAR*)rPtr;
                arrays1.k1Ptr = &rho;
                // wPtr = rDPtr * rPtr
                // rho += wPtr * rPtr
                gSum_host(&arrays1, &slave_userFunc_digPrecondSum);
            };

        #if(mpiTest)
            MPI_Request req1, req2;
            label flag1, flag2;
            scalar rd1, wApA1, delta1;

            rd1 = 0.0;
            MPI_Iallreduce(&rho, &rd1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req1);
            MPI_Test(&req1, &flag1, MPI_STATUSES_IGNORE);
        #endif


            if (nIter == 0)
            {
                // record first beta = 0
                betas.push_back(0.0);

                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        pPtr[cell]  = wPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = (SCALAR*)pPtr;
                    arrays1.A2Ptr = (SCALAR*)wPtr;
                    // pPtr = wPtr
                    vectorCopy_host(&arrays1);
                }

                // matrix_.Amul(niu, p, interfaceBouCoeffs_, interfaces_, cmpt);
                SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                        niu,
                                        p,
                                        interfaceBouCoeffs_,
                                        interfaces_,
                                        cmpt,
                                        _amul_parameter,
                                        _matrix_translate,
                                        _coarseLevel);

                wApA = 0.0;
                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        wApA +=  niuPtr[cell] * pPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A2Ptr = (SCALAR*)niuPtr;
                    arrays1.A3Ptr = (SCALAR*)pPtr;
                    arrays1.k1Ptr = &wApA;
                    // wApA += niuPtr * pPtr
                    gSum_host(&arrays1, &slave_userFunc_sumProd);
                }

                #if(mpiTest)
                {
                    wApA1 = 0.0;
                    MPI_Iallreduce(&wApA, &wApA1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req2);

                    MPI_Test(&req1, &flag1, MPI_STATUSES_IGNORE);
                    MPI_Waitall(1, &req1, MPI_STATUSES_IGNORE);

                    MPI_Test(&req2, &flag2, MPI_STATUSES_IGNORE);
                    MPI_Waitall(1, &req2, MPI_STATUSES_IGNORE);
                    rho = rd1;
                    wApA = wApA1;
                }
                #else
                {
                    scalar rd1_[2], rd2[2];
                    rd1_[0] = rho;
                    rd1_[1] = wApA;
                    if(UPstream::parRun())
                    {
                        MPI_Allreduce(&rd1_, &rd2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        rho   = rd2[0];
                        wApA  = rd2[1];
                    }
                }
                #endif
            }
            else
            {
                // matrix_.Amul(s, w, interfaceBouCoeffs_, interfaces_, cmpt);
                SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                        s,
                                        w,
                                        interfaceBouCoeffs_,
                                        interfaces_,
                                        cmpt,
                                        _amul_parameter,
                                        _matrix_translate,
                                        _coarseLevel);

                scalar delta = 0.0;
                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        delta +=  wPtr[cell] * sPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A2Ptr = (SCALAR*)wPtr;
                    arrays1.A3Ptr = (SCALAR*)sPtr;
                    arrays1.k1Ptr = &delta;
                    // delta += wPtr * sPtr
                    gSum_host(&arrays1, &slave_userFunc_sumProd);
                }

                #if(mpiTest)
                {
                    delta1 = 0.0;
                    MPI_Iallreduce(&delta, &delta1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req2);

                    MPI_Test(&req1, &flag1, MPI_STATUSES_IGNORE);
                    MPI_Waitall(1, &req1, MPI_STATUSES_IGNORE);
                    rho = rd1;
                }
                #else
                {
                    scalar rd1_[2], rd2[2];
                    rd1_[0] = rho;
                    rd1_[1] = delta;
                    if(UPstream::parRun())
                    {
                        MPI_Allreduce(&rd1_, &rd2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        rho   = rd2[0];
                        delta = rd2[1];
                    }
                }
                #endif




                scalar beta = rho / rhoOld;

                // record beta
                betas.push_back(beta);

                if(nCells < accUsingSize)
                {
                    for (register label cell=0; cell<nCells; cell++)
                    {
                        pPtr[cell]   = wPtr[cell] + beta * pPtr[cell];
                        niuPtr[cell] = sPtr[cell] + beta * niuPtr[cell];
                    }
                }
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr    = (SCALAR*)pPtr;
                    arrays1.A2Ptr    = (SCALAR*)niuPtr;
                    arrays1.A3Ptr    = (SCALAR*)wPtr;
                    arrays1.A4Ptr    = (SCALAR*)sPtr;
                    arrays1.k1       = beta;
                    arrays1.returnA2 = 1;
                    // pPtr   = wPtr + k1*pPtr
                    // niuPtr = sPtr + k1*niuPtr
                    vectorOps_host(&arrays1, &slave_userFunc_aEcPk1Mua_bEdPk1Mub);
                }

                #if(mpiTest)
                {
                    MPI_Test(&req2, &flag2, MPI_STATUSES_IGNORE);
                    MPI_Waitall(1, &req2, MPI_STATUSES_IGNORE);
                    delta = delta1;
                }
                #endif


                wApA  = delta - beta * beta * wApA;
            }

            // Update solution and residual:
            scalar alpha = rho / wApA;

            // record alpha
            alphas.push_back(alpha);

            if(nCells < accUsingSize)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    xPtr[cell] += alpha * pPtr[cell];
                    rPtr[cell] -= alpha * niuPtr[cell];
                }
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr    = (SCALAR*)xPtr;
                arrays1.A2Ptr    = (SCALAR*)rPtr;
                arrays1.A3Ptr    = (SCALAR*)pPtr;
                arrays1.A4Ptr    = (SCALAR*)niuPtr;
                arrays1.k1       = alpha;
                arrays1.returnA2 = 1;
                // xPtr += k1 * pPtr
                // rPtr -= k1 * niuPtr
                vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Muc_bEbMik1Mud);
            }
        }
        computeMaxEig(alphas,betas,nPCGsCheby_,1);
    }
    else
    {
        scalar maxEig, minEig;
        maxEig = maxEigPCG_ * boostFactorCheby_;
        minEig = maxEigPCG_ / eigRatioCheby_;

        scalar delta, theta, s1;
        scalar rhok, rhokp1, dtemp1, dtemp2;
        delta = 2.0 / (maxEig - minEig);
        theta = 0.5 * (maxEig + minEig);
        s1    = theta * delta;
        rhok  = 1.0 / s1;

        // Calculate A.x
        // matrix_.Amul(w, x, interfaceBouCoeffs_, interfaces_, cmpt);
        SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                w,
                                x,
                                interfaceBouCoeffs_,
                                interfaces_,
                                cmpt,
                                _amul_parameter,
                                _matrix_translate,
                                _coarseLevel);

        // At sweep = 0, update x
        scalar Rtheta = 1.0 / theta;
        if(nCells < accUsingSize)
        {
            for (register label cell=0; cell<nCells; cell++)
            {
                wPtr[cell]  = rDPtr[cell] * (bPtr[cell] - wPtr[cell]) * Rtheta;
                xPtr[cell] += wPtr[cell];
            }
        }
        else
        {
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr    = (SCALAR*)wPtr;
            arrays1.A2Ptr    = (SCALAR*)xPtr;
            arrays1.A3Ptr    = (SCALAR*)rDPtr;
            arrays1.A4Ptr    = (SCALAR*)bPtr;
            arrays1.k1       = Rtheta;
            arrays1.returnA2 = 1;
            // wPtr  = rDPtr * (bPtr - wPtr) * Rtheta
            // xPtr += wPtr
            vectorOps_host(&arrays1, &slave_userFunc_aEcMuSdMiaSMuk1_bEbPa);
        }

        for (label sweep=1; sweep<nSweeps; sweep++)
        {
            rhokp1 = 1.0 / (2.0 * s1 - rhok);
            dtemp1 = rhokp1 * rhok;
            dtemp2 = 2.0 * rhokp1 * delta;
            rhok   = rhokp1;

            // matrix_.Amul(p, x, interfaceBouCoeffs_, interfaces_, cmpt);
            SW_AMUL_COARSE_RELATION_REDESIGNED(matrix_,
                                    p,
                                    x,
                                    interfaceBouCoeffs_,
                                    interfaces_,
                                    cmpt,
                                    _amul_parameter,
                                    _matrix_translate,
                                    _coarseLevel);

            if(nCells < accUsingSize)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    wPtr[cell]  = dtemp1 * wPtr[cell] + dtemp2 * rDPtr[cell] * (bPtr[cell] - pPtr[cell]);
                    xPtr[cell] += wPtr[cell];
                }
            }
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (SCALAR*)wPtr;
                arrays1.A2Ptr = (SCALAR*)rDPtr;
                arrays1.A3Ptr = (SCALAR*)bPtr;
                arrays1.A4Ptr = (SCALAR*)pPtr;
                arrays1.k1    = dtemp1;
                arrays1.k2    = dtemp2;
                // wPtr = dtemp1 * wPtr + dtemp2 * rDPtr * (bPtr - pPtr)
                vectorOps_host(&arrays1, &slave_userFunc_aEk1MuaPk2MubMuScMidS);

                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = (SCALAR*)xPtr;
                arrays1.A2Ptr = (SCALAR*)wPtr;
                // xPtr += wPtr
                vectorOps_host(&arrays1, &slave_userFunc_aEaPb);
            }
        }

    }
}

void Foam::SW_ChebyIfpack2srPCG::smooth
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        fieldName_,
        x,
        matrix_,
        b,
        interfaceBouCoeffs_,
        interfaceIntCoeffs_,
        interfaces_,
        cmpt,
        nSweeps
    );
}
