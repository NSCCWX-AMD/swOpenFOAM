#include "GAMGAgglomeration.H"
#include "lduMatrix.H"
#include "labelField.H"
#include "primitiveFields.H"
#include "LUscalarMatrix.H"
#include "GAMGSolver.H"
#include "swGAMGSolver.hpp"
#include <mpi.h>
#define TIMES_SWAmul
#ifdef TIMES_SWAmul
#include "swAmulTranslate.hpp"
#include "amulMacros.h"
#endif
#include "swAmul.h"
extern "C" {
    #include "swAmul_host.h"
    #include "util.h"
}

namespace Foam
{

void SW_Amul
(
    const lduMatrix & matrix,
    scalarField& Apsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    amul_para &amul_parameter,
    amul_translate_array *matrix_translate
)
{
    scalar* __restrict__ ApsiPtr = Apsi.begin();
    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();

    // Initialise the update of interfaced interfaces
    matrix.initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

// std::cout << "Here Amul ok 0" << std::endl;
    amul_parameter.psiPtr = (SCALAR*)(&psiPtr[0]);
    amul_parameter.ApsiPtr =(SCALAR*)(&ApsiPtr[0]);
    ///init matrix_translate psi Apsi , by  lcx
    for(int i = 0 ; i < CORE_SIZE ; i++) {
        matrix_translate[i].psiPtr = amul_parameter.psiPtr;
        matrix_translate[i].ApsiPtr = amul_parameter.ApsiPtr;
    }

    amul_host(&amul_parameter, matrix_translate);

// std::cout << "Here Amul ok 1" << std::endl;

    // Update interface interfaces
    matrix.updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    tpsi.clear();
}

void SW_Tmul
(
    const lduMatrix & matrix,
    scalarField& Tpsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    amul_para &amul_parameter ,
    amul_translate_array *matrix_translate
)
{
    scalar* __restrict__ TpsiPtr = Tpsi.begin();

    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();


    // Initialise the update of interfaced interfaces
    matrix.initMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    amul_parameter.psiPtr = (SCALAR*)(&psiPtr[0]);
    amul_parameter.ApsiPtr =(SCALAR*)(&TpsiPtr[0]);
    ///init matrix_translate psi Apsi , by  lcx
    for(int i = 0 ; i < CORE_SIZE ; i++) {
        matrix_translate[i].psiPtr = amul_parameter.psiPtr;
        matrix_translate[i].ApsiPtr = amul_parameter.ApsiPtr;
    }
    tmul_host(&amul_parameter, matrix_translate);

    // Update interface interfaces
    matrix.updateMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    tpsi.clear();
}


void SW_ATmul
(
    const lduMatrix& matrix,
    scalarField& Apsi,
    scalarField& Tpsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    amul_para &amul_parameter,
    amul_translate_array *matrix_translate
)
{
    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();
    amul_parameter.psiPtr = (SCALAR*)(&psiPtr[0]);
    // Initialise the update of interfaced interfaces
    matrix.initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    const scalar* __restrict__ ApsiPtr = Apsi.begin();
    amul_parameter.ApsiPtr =(SCALAR*)(&ApsiPtr[0]);

    for(int i = 0 ; i < CORE_SIZE ; i++)
    {
        matrix_translate[i].psiPtr = amul_parameter.psiPtr;
        matrix_translate[i].ApsiPtr = amul_parameter.ApsiPtr;
    }
    amul_host(&amul_parameter, matrix_translate);

    const scalar* __restrict__ TpsiPtr = Tpsi.begin();
    amul_parameter.ApsiPtr = (SCALAR*)(&TpsiPtr[0]);

    for(int i = 0 ; i < CORE_SIZE ; i++)
    {
        matrix_translate[i].psiPtr = amul_parameter.psiPtr;
        matrix_translate[i].ApsiPtr = amul_parameter.ApsiPtr;
    }
    tmul_host(&amul_parameter, matrix_translate);

    // Update interface interfaces
    matrix.updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    matrix.updateMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );
    tpsi.clear();
}
} // end namespace Foam

