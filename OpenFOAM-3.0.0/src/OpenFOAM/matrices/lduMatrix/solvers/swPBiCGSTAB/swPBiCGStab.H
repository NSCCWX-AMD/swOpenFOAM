/*---------------------------------------------------------------------------*\
Class
    Foam::SW_PBiCGSTAB

Description
    Preconditioned bi-conjugate gradient stabilized solver for asymmetric
    lduMatrices using a run-time selectable preconditioner.

    References:
    \verbatim
        Van der Vorst, H. A. (1992).
        Bi-CGSTAB: A fast and smoothly converging variant of Bi-CG
        for the solution of nonsymmetric linear systems.
        SIAM Journal on scientific and Statistical Computing, 13(2), 631-644.

        Barrett, R., Berry, M. W., Chan, T. F., Demmel, J., Donato, J.,
        Dongarra, J., Eijkhout, V., Pozo, R., Romine, C. & Van der Vorst, H.
        (1994).
        Templates for the solution of linear systems:
        building blocks for iterative methods
        (Vol. 43). Siam.
    \endverbatim

SourceFiles
    SW_PBiCGSTAB.cpp

\*---------------------------------------------------------------------------*/

#ifndef SW_PBiCGStab_H
#define SW_PBiCGStab_H

#include "lduMatrix.H"
#include "swAmul.h"
#include "swAmulTranslate.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SW_PBiCGSTAB Declaration
\*---------------------------------------------------------------------------*/

class SW_PBiCGSTAB
:
    public lduMatrix::solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        SW_PBiCGSTAB(const SW_PBiCGSTAB&);

        //- Disallow default bitwise assignment
        void operator=(const SW_PBiCGSTAB&);

        const static int MAX_SW_USING_CORASE_LEVELS = 1;
    protected:
        static int if_first;
        static refilltion *_refill;

    public:

        static amul_translate_array **_matrix_traslate;

        static amul_para *_amul_parameter;
        static int _coarseLevel;


public:

    //- Runtime type information
    TypeName("swPBiCGStab");


    // Constructors

        //- Construct from matrix components and solver data stream
        SW_PBiCGSTAB
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~SW_PBiCGSTAB();
    // {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
