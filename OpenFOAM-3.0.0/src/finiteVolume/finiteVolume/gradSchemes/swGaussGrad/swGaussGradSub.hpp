#ifndef swGaussGradSub_H
#define swGaussGradSub_H

#include "vectorField.H"
#include "volMesh.H"
#include "volFieldsFwd.H"
#include "SlicedGeometricField.H"

namespace Foam
{
// the kenel template navigating to corresponding C routines
// GradType is the tensor product type of vector(Vector<scalar>)
// and Type. i.e.
// the possible match is
//      Type:       GradType:
//      scalar      vector
//      vector      tensor(Tensor<scalar>)
//      ... ...
// For detail, see class outerProduct and typeOfRank in Caelus.

/* declaration */
template<class Type, class GradType>
void swGaussGradInnerIteration
(const vectorField& Sf, const Field<Type>& issf, Field<GradType>& igGrad,
 const labelUList& owner, const labelUList& neighbour);

template<class GradType>
void swGaussGradDivide
(Field<GradType>& igGrad,
 const volScalarField::DimensionedInternalField& volume);

/* specialization declarations */
template<>
void swGaussGradInnerIteration<scalar, vector>
(const vectorField& Sf, const Field<scalar>& issf, Field<vector>& igGrad,
 const labelUList& owner, const labelUList& neighbour);

template<>
void swGaussGradInnerIteration<vector, tensor>
(const vectorField& Sf, const Field<vector>& issf, Field<tensor>& igGrad,
 const labelUList& owner, const labelUList& neighbour);

template<>
void swGaussGradDivide<vector>
(Field<vector>& igGrad,
 const volScalarField::DimensionedInternalField& volume);

template<>
void swGaussGradDivide<tensor>
(Field<tensor>& igGrad,
 const volScalarField::DimensionedInternalField& volume);

/* dummy implementation */
template<class Type, class GradType>
void swGaussGradInnerIteration
(const vectorField& Sf, const Field<Type>& issf, Field<GradType>& igGrad,
 const labelUList& owner, const labelUList& neighbour)
{
    Info<<"***Error:This template was not implemented!!!"<<endl;
    ::exit(1);
}

template<class GradType>
void swGaussGradDivide
(Field<GradType>& igGrad,
 const volScalarField::DimensionedInternalField& volume)
{
    Info<<"***Error:This template was not implemented!!!"<<endl;
    ::exit(1);
}

} // namespace Foam


#endif
