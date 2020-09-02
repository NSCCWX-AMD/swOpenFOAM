/*---------------------------------------------------------------------------*\
Class
    Foam::ChebyshevSmootherIfpack2

Description
    Chebyshev smoother with eigenvalue from PCG routine.

  Characteristics:
      - eigenvalue: using self-defined PCG routine
      - self-defined PCG routine: different from Caelus PCG,
        unessential parts removed
      -

SourceFiles
    ChebyshevSmootherIfpack2.cpp

Reference
    [1] Templates for the Solution of Linear Systems: Building Blocks
        for Iterative Methods, R. Barrett, M. Barry, T.F. Chan, J. Demmel,
        J. Donato, J. Dongarra, V. Eijkhout, R. Pozo, C. Romine, and
        Van der Vorst, SIAM, 1994, Philadephia, PA, 2nd edition

    [2] Left-scaled Chebyshev iteration, Ifpack2: Tempated Object-Oriented
        Algebraic Preconditioner Package, Mark Hoemmen, Sandia Corporation

    [3] On the use of conjugate gradient to calculate the eigenvalues and
        singular values of large, sparse matrices, John A. Scales, Geophysical
        Journal (1989) 97, 179-183

\*---------------------------------------------------------------------------*/


#ifndef ChebyshevSmootherIfpack2_H
#define ChebyshevSmootherIfpack2_H

#include "lduMatrix.H"
#include <vector>

namespace Foam
{
class ChebyshevSmootherIfpack2
:
    public lduMatrix::smoother
{
protected:
    //- Number of PCGs to get alphas and betas, used to calculate eigenvalue
    label nPCGsCheby_;

    //- Upper bound of the bounding ellipse of the eigenvalues of the matrix A
    mutable scalar maxEigPCG_;

    //- Detect if eigenvalue has been calculated,
    //- if not, Chebyshev smoother will call PCGs
    mutable bool eigFirstTimeComputed_;

    //- Read control parameters of Chebyshev smoother from dictionary
    void readDictCheby();

    //- Estimate the minEig by maxEig divided by eigRatio
    //- minEig = maxEig / eigRatioCheby_
    scalar eigRatioCheby_;

    //- factor to enlarge the maxEig, because
    //- the maxEig obtained is usually underestimated
    //- maxEig *= boostFactorCheby_
    scalar boostFactorCheby_;

    //- not used here
    scalar eigRatioCoarest_;

    //-detect if using preSmooth_
    //-if using, first Amul operation can be ignored
    label preSmoothUsing_;

    //- Sturm Sequence Method
    //- Find polynomial vector p(i) to find eigen value
    void h14Sturm
    (
        scalar **TriMatrix,         /// tridiagonal matrix
        const scalar lamb,          /// guess eigen value (used in bisection method)
        std::vector<scalar> &p,     /// a pvector polynomials (pi_1 to pi_n)
        label &s,                   /// the number of sign change
        const label nPCGs           /// nPCGsCheby_
    ) const;

    //- Form the matrix from alphas and betas, used to calculate eigenvalue
    //- algorithm see Ref.[3]
    void computeValueForMatrix
    (
        const std::vector<scalar> alphas,     /// alphas from nPCGs times of PCG
        const std::vector<scalar> betas,      /// betas from nPCGs times of PCG
        scalar **TriMatrix,              /// tridiagonal matrix
        const label nPCGs                /// nPCGsCheby_
    ) const;

    //- Estimate the range of eigenvalues
    void determineEigRange
    (
        scalar **TriMatrix,      /// tridiagonal matrix
        scalar &xBegin,          /// left bound of eigenvalues
        scalar &xEnd,            /// right bound of eigenvalues
        const label nPCGs        /// nPCGsCheby_
    ) const;

    //- Compute eigenvalue
    void computeMaxEig
    (
        const std::vector<scalar> alphas, /// alphas from nPCGs times of PCG
        const std::vector<scalar> betas,  /// betas from nPCGs times of PCG
        const label nPCGs,           /// nPCGsCheby_
        const label k                /// kth largest eigenvalue, default k=1 for the maximum eigenvalue
    ) const;

    template <typename T> T **allocateSym2D(label nPCGs) const;
    template <typename T> void deleteSym2D(T **arr, label nPCGs) const;

public:

    //- Runtime type information
    TypeName("Cheby");

    // Constructors

        //- Construct from components
        ChebyshevSmootherIfpack2
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Destructor
        virtual ~ChebyshevSmootherIfpack2();


    // Member Functions

        //- Smooth for the given number of sweeps
        // static void smooth
        void smooth
        (
            const word& fieldName,
            scalarField& x,
            const lduMatrix& matrix,
            const scalarField& b,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt,
            const label nSweeps
        ) const;


        //- Smooth the solution for a given number of sweeps
        virtual void smooth
        (
            scalarField& x,
            const scalarField& B,
            const direction cmpt,
            const label nSweeps
        ) const;
};


} // End namespace Foam

#endif
