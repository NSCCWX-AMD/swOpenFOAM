#ifndef SW_ChebyIfpack2srPCG_H
#define SW_ChebyIfpack2srPCG_H

#include "lduMatrix.H"
#include <vector>
#include "ChebyshevSmootherIfpack2.hpp"
#include "amulMacros.h"
namespace Foam
{

class SW_ChebyIfpack2srPCG
    :public ChebyshevSmootherIfpack2
{

private:
    int _coarseLevel;
    amul_para * _amul_parameter ;
    amul_translate_array * _matrix_translate;
    static int MAX_SW_USING_CORASE_LEVELS;
public:

    //- Runtime type information
    TypeName("swCheby");

    // Constructors
        static void  SET_MAX_SW_USING_CORASE_LEVELS(int level);
        //- Construct from components
        SW_ChebyIfpack2srPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Destructor
        virtual ~SW_ChebyIfpack2srPCG();


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
