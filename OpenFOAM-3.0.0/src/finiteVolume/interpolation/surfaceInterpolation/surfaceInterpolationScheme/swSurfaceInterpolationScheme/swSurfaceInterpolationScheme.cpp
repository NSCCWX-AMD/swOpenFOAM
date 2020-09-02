/*---------------------------------------------------------------------------*\
Copyright (C) 2011 OpenFOAM Foundation
-------------------------------------------------------------------------------
surfaceInterpolationScheme.cpp contains specializtions of
interpolate members of surfaceInterpolationScheme calss
in surfaceInterpolationScheme.hpp to use Sunway slave cores

author: Hu Ren
email: renhu@mail.nsccwx.cn
\*---------------------------------------------------------------------------*/

#include "surfaceInterpolationScheme.H"
#include "swSurfaceInterpolationScheme_C.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// 1st overload:
//- Return the face-interpolate of the given cell field
//  with the given owner and neighbour weighting factors
// 2nd overload:
//- Return the face-interpolate of the given cell field
//  with the given weigting factors
#define specializeInterpolate(Type) \
template<>\
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >\
surfaceInterpolationScheme<Type >::interpolate\
(\
    const GeometricField<Type, fvPatchField, volMesh>& vf,\
    const tmp<surfaceScalarField>& tlambdas,\
    const tmp<surfaceScalarField>& tys\
)\
{\
    if (surfaceInterpolation::debug)\
    {\
        Info<< "surfaceInterpolationScheme<Type >::uncorrectedInterpolate"\
               "(const GeometricField<Type, fvPatchField, volMesh>&, "\
               "const tmp<surfaceScalarField>&, "\
               "const tmp<surfaceScalarField>&) : "\
               "interpolating volTypeField from cells to faces "\
               "without explicit correction"\
            << endl;\
    }\
\
    const surfaceScalarField& lambdas = tlambdas();\
    const surfaceScalarField& ys = tys();\
\
    const Field<Type >& vfi = vf.internalField();\
    const scalarField& lambda = lambdas.internalField();\
    const scalarField& y = ys.internalField();\
\
    const fvMesh& mesh = vf.mesh();\
    const labelUList& P = mesh.owner();\
    const labelUList& N = mesh.neighbour();\
\
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf\
    (\
        new GeometricField<Type, fvsPatchField, surfaceMesh>\
        (\
            IOobject\
            (\
                "interpolate("+vf.name()+')',\
                vf.instance(),\
                vf.db()\
            ),\
            mesh,\
            vf.dimensions()\
        )\
    );\
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();\
\
    Field<Type >& sfi = sf.internalField();\
\
    interpolate1(\
        P.size(),\
        vfi.size(),\
        &sfi[0],        \
        &vfi[0],         \
        &lambda[0],          \
        &y[0],           \
        &P[0],           \
        &N[0],          \
        sizeof(sfi[0])\
        ); \
\
    forAll(lambdas.boundaryField(), pi)\
    {\
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];\
        const fvsPatchScalarField& pY = ys.boundaryField()[pi];\
\
        if (vf.boundaryField()[pi].coupled())\
        {\
            sf.boundaryField()[pi] =\
                pLambda*vf.boundaryField()[pi].patchInternalField()\
              + pY*vf.boundaryField()[pi].patchNeighbourField();\
        }\
        else\
        {\
            sf.boundaryField()[pi] = vf.boundaryField()[pi];\
        }\
    }\
\
    tlambdas.clear();\
    tys.clear();\
\
    return tsf;\
} \
 \
 \
template<>\
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >\
surfaceInterpolationScheme<Type >::interpolate\
(\
    const GeometricField<Type, fvPatchField, volMesh>& vf,\
    const tmp<surfaceScalarField>& tlambdas\
)\
{\
    if (surfaceInterpolation::debug)\
    {\
        Info<< "surfaceInterpolationScheme<Type >::interpolate"\
               "(const GeometricField<Type, fvPatchField, volMesh>&, "\
               "const tmp<surfaceScalarField>&) : "\
               "interpolating volTypeField from cells to faces "\
               "without explicit correction"\
            << endl;\
    }\
\
    const surfaceScalarField& lambdas = tlambdas();\
\
    const Field<Type >& vfi = vf.internalField();\
    const scalarField& lambda = lambdas.internalField();\
\
    const fvMesh& mesh = vf.mesh();\
    const labelUList& P = mesh.owner();\
    const labelUList& N = mesh.neighbour();\
\
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf\
    (\
        new GeometricField<Type, fvsPatchField, surfaceMesh>\
        (\
            IOobject\
            (\
                "interpolate("+vf.name()+')',\
                vf.instance(),\
                vf.db()\
            ),\
            mesh,\
            vf.dimensions()\
        )\
    );\
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();\
\
    Field<Type >& sfi = sf.internalField();\
\
    interpolate1(\
        P.size(),\
        vfi.size(),\
        sfi.cdata(),        \
        vfi.cdata(),         \
        lambda.cdata(),          \
        lambda.cdata(),           \
        P.cdata(),           \
        N.cdata(),          \
        sizeof(sfi[0])\
        );\
\
    forAll(lambdas.boundaryField(), pi)\
    {\
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];\
\
        if (vf.boundaryField()[pi].coupled())\
        {\
            tsf().boundaryField()[pi] =\
                pLambda*vf.boundaryField()[pi].patchInternalField()\
             + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();\
        }\
        else\
        {\
            sf.boundaryField()[pi] = vf.boundaryField()[pi];\
        }\
    }\
\
    tlambdas.clear();\
\
    return tsf;\
}

#ifdef SW_SLAVE
specializeInterpolate(double)
specializeInterpolate(Vector<double>)
specializeInterpolate(SphericalTensor<double>)
specializeInterpolate(SymmTensor<double>)
specializeInterpolate(Tensor<double>)
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
