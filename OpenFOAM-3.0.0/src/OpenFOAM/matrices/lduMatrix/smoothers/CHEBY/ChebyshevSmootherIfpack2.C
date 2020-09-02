#include "ChebyshevSmootherIfpack2.hpp"
#include "DICPreconditioner.H"
#include "smoothSolver.H"
#include <time.h>
#include "diagonalPreconditioner.H"
#include <vector>


namespace Foam
{
    defineTypeNameAndDebug(ChebyshevSmootherIfpack2, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<ChebyshevSmootherIfpack2>
        addChebyshevSmootherIfpack2SymMatrixConstructorToTable_;
}

Foam::ChebyshevSmootherIfpack2::ChebyshevSmootherIfpack2
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    nPCGsCheby_(10),
    maxEigPCG_(0.0),
    eigFirstTimeComputed_(true),
    eigRatioCheby_(1.5),
    boostFactorCheby_(1.1),
    eigRatioCoarest_(1.0),
    preSmoothUsing_(0)
{
    readDictCheby();
}

Foam::ChebyshevSmootherIfpack2::~ChebyshevSmootherIfpack2()
{
}

void Foam::ChebyshevSmootherIfpack2::readDictCheby()
{
    // we could also consider supplying defaults here too
    controlDict_.readIfPresent("nPCGsCheby",       nPCGsCheby_);
    controlDict_.readIfPresent("eigRatioCheby",    eigRatioCheby_);
    controlDict_.readIfPresent("boostFactorCheby", boostFactorCheby_);
    controlDict_.readIfPresent("eigRatioCoarest",  eigRatioCoarest_);
    controlDict_.readIfPresent("nPreSweeps",       preSmoothUsing_);
}

void Foam::ChebyshevSmootherIfpack2::h14Sturm
(
    scalar **TriMatrix,
    const scalar lamb,
    std::vector<scalar> &p,
    label &s,
    const label nPCGs
) const
{
    label k;
    s = 0;

    p[0] = lamb - TriMatrix[0][0];
    if (p[0] < 0)
        s = 1;

    p[1] = (lamb - TriMatrix[1][1]) * p[0] -pow(TriMatrix[0][1], 2);
    if(p[1] * p[0] < 0)
        s++;

    for(k=2; k<nPCGs; k++)
    {
        p[k] = (lamb - TriMatrix[k][k]) * p[k-1] - pow(TriMatrix[k-1][k], 2)* p[k-2];
        if(p[k] * p[k-1] < 0)
            s++;
    }
}

// allocate memory
template <typename T>
T** Foam::ChebyshevSmootherIfpack2::allocateSym2D (label n) const
{
    T **arr2D;
    arr2D = new T *[n];

    for(label i=0; i<n; i++)
    {
        arr2D[i] = new T[n];
    }

    return (T**)arr2D;
}


// free memory
template <typename T>
void Foam::ChebyshevSmootherIfpack2::deleteSym2D(T **arr, label n) const
{
	if(arr != NULL)
	{
        for(label i=0; i<n; i++)
		{
		    delete [] arr[i];
			arr[i] = NULL;
		}

		delete [] arr;
        arr = NULL;
	}      
}


void Foam::ChebyshevSmootherIfpack2::computeValueForMatrix
(
    const std::vector<scalar> alphas,
    const std::vector<scalar> betas,
    scalar **TriMatrix,
    const label nPCGs
) const
{
    for (label i=0; i<nPCGs; i++)
        for (label j=0; j<nPCGs; j++)
        {
            TriMatrix[i][j] = 0;
        }

    for (label i=0; i<nPCGs; i++)
    {
        // diagonal values
        if(i == 0)
            TriMatrix[i][i] = 1.0/alphas[i];
        else
            TriMatrix[i][i] = 1.0/alphas[i] + betas[i]/alphas[i-1];

        // off-diag values
        if(i < nPCGs-1)
            TriMatrix[i][i+1] = sqrt(betas[i+1])/alphas[i];
        if(i > 0)
            TriMatrix[i][i-1] = TriMatrix[i-1][i];
    }
}

void Foam::ChebyshevSmootherIfpack2::determineEigRange
(
    scalar **TriMatrix,
    scalar &xBegin,
    scalar &xEnd,
    const label nPCGs
) const
{
    scalar Tright, Tleft;
    scalar lambMax, lambMin;

    lambMin = lambMax = xBegin = xEnd = 0.0;
    for(label i=0; i<nPCGs; i++)
    {
        Tright = (i == (nPCGs-1))? 0 : TriMatrix[i][i+1];
        Tleft  = (i == 0)?         0 : TriMatrix[i][i-1];
        lambMax = TriMatrix[i][i] + fabs(Tright) + fabs(Tleft);
        lambMin = TriMatrix[i][i] - fabs(Tright) - fabs(Tleft);
        xBegin = (xBegin >= lambMin)? lambMin : xBegin;
        xEnd   = (xEnd   <= lambMax)? lambMax : xEnd;
    }

}

void Foam::ChebyshevSmootherIfpack2::computeMaxEig
(
    const std::vector<scalar> alphas,
    const std::vector<scalar> betas,
    const label nPCGs,
    const label k
) const
{
    scalar xBegin, xEnd, xn;
    label  s;
    const scalar SMALL = 0.00001;
    std::vector<scalar> pvector(nPCGs);

    scalar **TriMatrix = allocateSym2D<scalar>(nPCGs);
    computeValueForMatrix(alphas, betas, TriMatrix, nPCGs);
    determineEigRange(TriMatrix, xBegin, xEnd, nPCGs);

    do
    {
        xn = (xBegin + xEnd) * 0.5;
        h14Sturm(TriMatrix, xn, pvector, s, nPCGs);
        pvector.insert(pvector.begin(), 1, 1);   // because p[0] = 1

        if(s >= k)
            xBegin = xn;
        else
            xEnd = xn;
    }while
    (
        (xEnd - xBegin) >= SMALL
    );

    deleteSym2D(TriMatrix, nPCGs);

    maxEigPCG_ = xn;
}

void Foam::ChebyshevSmootherIfpack2::smooth
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

    // diagonal Preconditioner
    diagonalPreconditioner diagPre
    (
        smoothSolver
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            dictionary()
        ),
        dictionary()
    );


    register const label nCells = x.size();
    scalar* __restrict__ xPtr = x.begin();

    scalarField p(nCells);
    scalar* __restrict__ pPtr = p.begin();

    scalarField w(nCells);
    scalar* __restrict__ wPtr = w.begin();

    scalarField r(nCells);
    scalar* __restrict__ rPtr = r.begin();

    const scalar* __restrict__ bPtr = b.begin();

    // compute eigenvalue using diagnoal conjugate gradient sweeps
    if(eigFirstTimeComputed_)
    {
        eigFirstTimeComputed_ = false;

        std::vector<scalar> alphas;
        std::vector<scalar> betas;

        scalar rho = 1e20;
        scalar rhoOld = rho;

        // if(preSmoothUsing_)
        // {
        //     for(label cell=0; cell<nCells; cell++)
        //     {
        //         rPtr[cell] = bPtr[cell];
        //     }
        // }
        // else
        {
            matrix_.Amul(w, x, interfaceBouCoeffs_, interfaces_, cmpt);

            for (label cell=0; cell<nCells; ++cell)
            {
                rPtr[cell] = bPtr[cell] - wPtr[cell];
            }
        }

        // Clear the alpha beta vectors
        while(alphas.begin() != alphas.end())
            alphas.pop_back();
        while( betas.begin() != betas.end() )
            betas.pop_back();

        for(label nIter=0; nIter<nPCGsCheby_; nIter++)
        {
            // Store previous rho
            rhoOld = rho;

            // Apply preconditioner
            diagPre.precondition(w, r, cmpt);

            // Update rho:
            rho = gSumProd(w, r);

            if (nIter == 0)
            {
                // Hanfeng
                // record first beta = 0
                betas.push_back(0.0);
                // Info << "At iterative step " << solverPerf.nIterations() << ", beta is " << 0 << endl;
                // Hanfeng

                for (register label cell=0; cell<nCells; cell++)
                {
                    pPtr[cell] = wPtr[cell];
                }
            }
            else
            {
                scalar beta = rho/rhoOld;

                // record beta
                betas.push_back(beta);

                for (register label cell=0; cell<nCells; cell++)
                {
                    pPtr[cell] = wPtr[cell] + beta*pPtr[cell];
                }
            }


            // Update preconditioned residual
            matrix_.Amul(w, p, interfaceBouCoeffs_, interfaces_, cmpt);

            scalar wApA = gSumProd(w, p);

            // Update solution and residual:
            scalar alpha = rho/wApA;

            // record alpha
            alphas.push_back(alpha);
            // Info << "At iterative step " << solverPerf.nIterations() << ", alpha is " << alpha << endl;

            for (register label cell=0; cell<nCells; cell++)
            {
                xPtr[cell] += alpha*pPtr[cell];
                rPtr[cell] -= alpha*wPtr[cell];
            }
        }

        computeMaxEig(alphas, betas, nPCGsCheby_, 1);
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
        matrix_.Amul(w, x, interfaceBouCoeffs_, interfaces_, cmpt);

        // Calculate initial residual field
        for (register label cell=0; cell<nCells; cell++)
        {
            rPtr[cell] = bPtr[cell] - wPtr[cell];
        }

        // At sweep = 0, update x
        diagPre.precondition(w, r, cmpt);

        for (register label cell=0; cell<nCells; cell++)
        {
            wPtr[cell]  = wPtr[cell] / theta;
            xPtr[cell] += wPtr[cell];
        }

        for (label sweep=1; sweep<nSweeps; sweep++)
        {
            rhokp1 = 1.0 / (2.0 * s1 - rhok);
            dtemp1 = rhokp1 * rhok;
            dtemp2 = 2.0 * rhokp1 * delta;
            rhok   = rhokp1;

            matrix_.Amul(p, x, interfaceBouCoeffs_, interfaces_, cmpt);

            for (register label cell=0; cell<nCells; cell++)
            {
                rPtr[cell] = bPtr[cell] - pPtr[cell];
            }

            diagPre.precondition(p, r, cmpt);

            for (register label cell=0; cell<nCells; cell++)
            {
                wPtr[cell]  = dtemp1 * wPtr[cell] + dtemp2 * pPtr[cell];
                xPtr[cell] += wPtr[cell];
            }
        }
    }
}

void Foam::ChebyshevSmootherIfpack2::smooth
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
