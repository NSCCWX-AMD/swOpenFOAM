/*---------------------------------------------------------------------------*\
Copyright (C) 2011 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of CAELUS.

    CAELUS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CAELUS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CAELUS.  If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam::fvc

Description
    Surface integrate surfaceField creating a volField.
    Surface sum a surfaceField creating a volField.

SourceFiles
    fvcSurfaceIntegrate.cpp

author：zhangyongbin
add SPE algorithm access:packageData in fvcSurfaceIntegrateSub.hpp


\*---------------------------------------------------------------------------*/


//#include "volFieldsFwd.hpp"
//#include "surfaceFieldsFwd.hpp"
//#include "fvMesh.hpp"
//#include "zeroGradientFvPatchFields.hpp"
#include "fvcSurfaceIntegrateSub.hpp"
#include "fvsPatchField.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "fvcSurfaceIntegrate.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// This is the <double> type specialization
template<>
void surfaceIntegrate<double>
(
    Field<double>& ivf,
    const GeometricField<double, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field<double>& issf = ssf;

#ifdef SW_SLAVE
 {

     packageData(owner,neighbour,issf,ivf);

 }
#else
    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
#endif
    //Info<<ivf.first()<<"<->"<<ivf.last()<<endl;
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<double>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();
}

// This is the <vector> type specialization
template<>
void surfaceIntegrate< Vector<double> >
(
    Field< Vector<double> >& ivf,
    const GeometricField< Vector<double>, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field< Vector<double> >& issf = ssf;

#ifdef SW_SLAVE
 {

     packageData(owner,neighbour,issf,ivf);

 }
#else
    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
#endif
    //Info<<ivf.first()<<"<->"<<ivf.last()<<endl;
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField< Vector<double> >& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();
}

// This is the <tensor> type specialization
template<>
void surfaceIntegrate< Tensor<double> >
(
    Field< Tensor<double> >& ivf,
    const GeometricField< Tensor<double>, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field< Tensor<double> >& issf = ssf;

#ifdef SW_SLAVE
 {

     packageData(owner,neighbour,issf,ivf);

 }
#else
    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
#endif
    //Info<<ivf.first()<<"<->"<<ivf.last()<<endl;
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField< Tensor<double> >& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();
}

// This is the <symmTensor> type specialization
template<>
void surfaceIntegrate< SymmTensor<double> >
(
    Field< SymmTensor<double> >& ivf,
    const GeometricField< SymmTensor<double>, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field< SymmTensor<double> >& issf = ssf;

#ifdef SW_SLAVE
 {

     packageData(owner,neighbour,issf,ivf);

 }
#else
    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
#endif
    //Info<<ivf.first()<<"<->"<<ivf.last()<<endl;
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField< SymmTensor<double> >& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();
}

// This is the <sphericalTensor> type specialization
template<>
void surfaceIntegrate< SphericalTensor<double> >
(
    Field< SphericalTensor<double> >& ivf,
    const GeometricField< SphericalTensor<double>, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field< SphericalTensor<double> >& issf = ssf;

#ifdef SW_SLAVE
 {

     packageData(owner,neighbour,issf,ivf);

 }
#else
    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
#endif
    //Info<<ivf.first()<<"<->"<<ivf.last()<<endl;
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField< SphericalTensor<double> >& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
