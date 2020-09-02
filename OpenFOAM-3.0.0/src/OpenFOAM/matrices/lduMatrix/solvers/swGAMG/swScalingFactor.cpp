#include "swGAMGSolver.hpp"
#include "vector2D.H"
#include "amulMacros.h"
#include "vectorOps.h"
#include "mpi.h"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::SW_GAMGSolver::scalingFactor
(
    scalarField& Acf,
    const lduMatrix& A,
    scalarField& field,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const scalarField& source,
    const direction cmpt,
    const int coarseLevel
) const
{
    // A.Amul(Acf, field, interfaceLevelBouCoeffs, interfaceLevel, cmpt);
    SW_AMUL_COARSE_RELATION_REDESIGNED(A,
                            Acf,
                            field,
                            interfaceLevelBouCoeffs,
                            interfaceLevel,
                            cmpt,
                            &(_amul_parameter[coarseLevel]),
                            _matrix_traslate[coarseLevel],
                            coarseLevel);
    const label nCells = field.size();
    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;
    const scalarField &D = A.diag();

    scalar vec_temp1[2], vec_temp2[2];
    if(nCells < 2500)
    {
        for (register label cell=0; cell<nCells; cell++)
        {
            scalingFactorNum   += source[cell]*field[cell];
            scalingFactorDenom += Acf   [cell]*field[cell];
            // While the matrix-multiply done for the scaling it is
            // possible to perform a point-Jacobi smoothing operation cheaply
            field[cell] += (source[cell] - Acf[cell]) / D[cell];
        }
        vec_temp1[0] = scalingFactorNum;
        vec_temp1[1] = scalingFactorDenom;
    }
    else
    {
        MVM_Arrays arrays1;
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = (SCALAR*)&(field[0]);
        arrays1.A2Ptr = (SCALAR*)&(source[0]);
        arrays1.A3Ptr = (SCALAR*)&(Acf[0]);
        arrays1.A4Ptr = (SCALAR*)&(D[0]);
        arrays1.k1Ptr = &vec_temp1[0];
        arrays1.k2Ptr = &vec_temp1[1];
        scalingFactor_host(&arrays1);
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

    return vec_temp2[0]/stabilise(vec_temp2[1], VSMALL);

    // vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    // reduce(scalingVector, sumOp<vector2D>());
    // return scalingVector.x()/stabilise(scalingVector.y(), VSMALL);
}


