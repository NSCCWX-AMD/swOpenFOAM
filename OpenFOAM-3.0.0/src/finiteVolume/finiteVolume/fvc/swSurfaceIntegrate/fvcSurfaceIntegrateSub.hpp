#ifndef fvcSurfaceIntegrateSub_H
#define fvcSurfaceIntegrateSub_H
#include "surfaceIntegrate_host.h"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fvc
{
    template<class Type>
    void packageData
    (
            const labelUList& ,const labelUList& ,const Field<Type>& , Field<Type>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template< >
void packageData<double>
(
     const labelUList& owner,
     const labelUList& neighbour,
     const Field<double>& issf,
     Field<double>& ivf
     )
{
     const int* ownerPtr=(const int*)owner.cdata();
     const int* neighbourPtr=(const int*)neighbour.cdata();
     const double* issfPtr=(const double*)issf.cdata();
     double* ivfPtr=(double *)ivf.cdata();

     const int issfSize=issf.size();
     const int ivfSize=ivf.size();
     const int owner_size=owner.size();
     const int neighbour_size=neighbour.size();
     const int vectorSize=1;

     surfaceIntegrate_host(ownerPtr,neighbourPtr,issfPtr,ivfPtr,issfSize,ivfSize,vectorSize,owner_size,neighbour_size);

}


template< >
void packageData< Vector<double> >
(
     const labelUList& owner,
     const labelUList& neighbour,
     const Field< Vector<double> >& issf,
     Field< Vector<double> >& ivf
     )
{
     const int* ownerPtr=owner.cdata();
     const int* neighbourPtr=neighbour.cdata();
     const  Vector<double>* issfPtr=issf.cdata();
     const  Vector<double>* ivfPtr=ivf.cdata();

     const double* issf_Ptr=(const double*)&issfPtr[0];
     double* ivf_Ptr=(double*)&ivfPtr[0];
     const int issfSize=issf.size();
     const int ivfSize=ivf.size();const int owner_size=owner.size();
     const int neighbour_size=neighbour.size();
     const int vectorSize=3;

     surfaceIntegrate_host(ownerPtr,neighbourPtr,issf_Ptr,ivf_Ptr,issfSize,ivfSize,vectorSize,owner_size,neighbour_size);
}


template< >
void packageData< Tensor<double> >
(
     const labelUList& owner,
     const labelUList& neighbour,
     const Field< Tensor<double> >& issf,
     Field< Tensor<double> >& ivf
     )
{
     const int* ownerPtr=owner.cdata();
     const int* neighbourPtr=neighbour.cdata();
     const  Tensor<double>* issfPtr=issf.cdata();
     const  Tensor<double>* ivfPtr=ivf.cdata();
     const double* issf_Ptr=(const double*)&issfPtr[0];
     double* ivf_Ptr=(double*)&ivfPtr[0];

     const int issfSize=issf.size();
     const int ivfSize=ivf.size();
     const int owner_size=owner.size();
     const int neighbour_size=neighbour.size();
     const int vectorSize=6;
    surfaceIntegrate_host(ownerPtr,neighbourPtr,issf_Ptr,ivf_Ptr,issfSize,ivfSize,vectorSize,owner_size,neighbour_size);
}


template< >
void packageData< SymmTensor<double> >
(
     const labelUList& owner,
     const labelUList& neighbour,
     const Field< SymmTensor<double> >& issf,
     Field< SymmTensor<double> >& ivf
     )
{
     const int* ownerPtr=owner.cdata();
     const int* neighbourPtr=neighbour.cdata();
     const  SymmTensor<double>* issfPtr=issf.cdata();
     const  SymmTensor<double>* ivfPtr=ivf.cdata();
     const double* issf_Ptr=(const double*)issfPtr;
     double* ivf_Ptr=(double *)ivfPtr;

     const int issfSize=issf.size();
     const int ivfSize=ivf.size();
     const int owner_size=owner.size();
     const int neighbour_size=neighbour.size();
     const int vectorSize=9;
     surfaceIntegrate_host(ownerPtr,neighbourPtr,issf_Ptr,ivf_Ptr,issfSize,ivfSize,vectorSize,owner_size,neighbour_size);
}


template< >
void packageData< SphericalTensor<double> >
(
     const labelUList& owner,
     const labelUList& neighbour,
     const Field< SphericalTensor<double> >& issf,
     Field< SphericalTensor<double> >& ivf
     )
{
     const int* ownerPtr=owner.cdata();
     const int* neighbourPtr=neighbour.cdata();
     const  SphericalTensor<double>* issfPtr=issf.cdata();
     const  SphericalTensor<double>* ivfPtr=ivf.cdata();
     const double* issf_Ptr=(const double*)issfPtr;
     double* ivf_Ptr=(double*)ivfPtr;

     const int issfSize=issf.size();
     const int ivfSize=ivf.size();
     const int owner_size=owner.size();
     const int neighbour_size=neighbour.size();
     const int vectorSize=1;
     surfaceIntegrate_host(ownerPtr,neighbourPtr,issf_Ptr,ivf_Ptr,issfSize,ivfSize,vectorSize,owner_size,neighbour_size);
}

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
