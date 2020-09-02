#ifndef SW_Amul_H
#define SW_Amul_H
#include "GAMGAgglomeration.H"
#include "lduMatrix.H"
#include "labelField.H"
#include "primitiveFields.H"
#include "LUscalarMatrix.H"
#include "GAMGSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "amulMacros.h"
namespace Foam{
	void SW_Amul
    (
        const lduMatrix & matrix,
        scalarField& Apsi,
        const tmp<scalarField>& tpsi,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const direction cmpt,
        amul_para &amul_parameter ,
        amul_translate_array *matrix_translate
    );
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
	);

    void SW_ATmul
    (
        const lduMatrix & matrix,
        scalarField& Apsi,
        scalarField& Tpsi,
        const tmp<scalarField>& tpsi,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const direction cmpt,
        amul_para &amul_parameter ,
        amul_translate_array *matrix_translate
    );

    void swAmulRedesigned
    (
        const lduMatrix & matrix,
        scalarField& Apsi,
        const tmp<scalarField>& tpsi,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const direction cmpt,
        amul_para &amul_parameter ,
        amul_translate_array *matrix_translate
    );

    void swGAMGAmulRedesigned
    (
        const lduMatrix & matrix,
        scalarField& Apsi,
        const tmp<scalarField>& tpsi,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const direction cmpt,
        amul_para &amul_parameter ,
        amul_translate_array *matrix_translate
    );

    void GAMGAmulRedesigned
    (
        const lduMatrix & matrix,
        scalarField& Apsi,
        const tmp<scalarField>& tpsi,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const direction cmpt
    );

    void initAndUpdate
    (
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const tmp<scalarField>& tpsi,
        scalarField& Apsi,
        const direction cmpt
    );

}
#endif
