#include "stdlib.h"
#include "rowSubsections.hpp"
#include "tensor.H"
#include "swGaussGradSub.hpp"
#include "swGgradfRoutines.h"

#include "mpi.h"

// macros used for serial debug in mpi parallel program
#define SERIAL_BEGIN(procNum) \
{ \
int mpi_id; \
MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id); \
for(int id = 0; id < procNum; id++) \
{ \
    if( id  == mpi_id) \
    { \
        printf("mpi rank %d:\n", id);

#define SERIAL_END \
    } \
    MPI_Barrier(MPI_COMM_WORLD); \
} \
}

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
template<>
void swGaussGradInnerIteration<scalar, vector>
    (const vectorField& Sf, const Field<scalar>& issf, Field<vector>& igGrad,
     const labelUList& owner, const labelUList& neighbour)
{

    const char* Sf_vptr    = (const char*) (Sf.cdata()); // char pointer pointing to a vector field
    static const label fSize      = Sf.size(); // face field size
    if(fSize != Sf.size()) Info<<"***Error: face array length changed from "
                               <<fSize<<" to "<< Sf.size()<<endl;
    const scalar* issf_ptr = issf.cdata();

    const label* lPtr = owner.cdata();
    const label* uPtr = neighbour.cdata();
    if( fSize != issf.size()  ||
        fSize != owner.size() ||
        fSize != neighbour.size())
    {
        Info<<endl<<"***Error: Face fields sizes are not the same!"<<endl<<endl;
        ::exit(1);
    }

    char* igGrad_vptr  = (char*) (igGrad.data());
    static const label cSize  = igGrad.size(); // cell field size
    if(cSize != igGrad.size()) Info<<"***Error: cell array length changed from "
                                   <<cSize<<" to "<<igGrad.size()<<endl;

    const label vector_size   = (label) sizeof(vector);
    const label vector_offset = (char*)(Sf[0].v_)-(char*)(&Sf[0]);


#if FULL_DEGUG
    Info<<endl<<"!! Init scalar RowSubsections for GaussGradf"<<endl<<endl;
#endif
    // calculate the max face number for LDM to include all relative arrays
    static label faceSize4CPE = (label)
                                (
                                    0.6*64.0*1000.0/
                                    (
										// data on edges
                                        (scalar)sizeof(scalar)
                                        + (scalar)sizeof(label)
                                        + (scalar)vector_size
										// data on columns, use edge length
                                        + (scalar)sizeof(scalar)*3.0
                                        + (scalar)vector_size
                                    )
                                );

    // static RowSubsections object that will not be construct repeatedly
    // when this function is called for more than once
    static SWFoam::RowSubsections secs(fSize, cSize, 64, 4, lPtr, uPtr,
        faceSize4CPE, vector_size);

    // call the C sub-routine
    swGgradfInnerRoutine_scalar( Sf_vptr, issf_ptr, igGrad_vptr, lPtr, uPtr,
                                 fSize, cSize, vector_size, vector_offset,
                                 secs.getSubsections(),secs.getSecNumInSeg(),
                                 secs.getColRoundNum());

}

template<>
void swGaussGradInnerIteration<vector, tensor>
    (const vectorField& Sf, const Field<vector>& issf, Field<tensor>& igGrad,
     const labelUList& owner, const labelUList& neighbour)
{
    const char* Sf_vptr    = (const char*) (Sf.cdata()); // char pointer pointing to a vector field
    static const label fSize      = Sf.size(); // face field size
    if(fSize != Sf.size()) Info<<"***Error: face array length changed from "
                               <<fSize<<" to "<< Sf.size()<<endl;
    const char* issf_vptr  = (const char*) (issf.cdata()); // char pointer pointing to a vector field

    const label* lPtr = owner.cdata();
    const label* uPtr = neighbour.cdata();
    if( fSize != issf.size()  ||
        fSize != owner.size() ||
        fSize != neighbour.size())
    {
        Info<<endl<<"***Error: Face fields sizes are not the same!"<<endl<<endl;
        ::exit(1);
    }

    char* igGrad_tptr  = (char*) (igGrad.data());
    static const label cSize  = igGrad.size(); // cell field size
    if(cSize != igGrad.size()) Info<<"***Error: cell array length changed from "
                                   <<cSize<<" to "<<igGrad.size()<<endl;

    const label vector_size   = (label) sizeof(vector);
    const label vector_offset = (char*)(Sf[0].v_)-(char*)(&Sf[0]);

    const label tensor_size   = (label) sizeof(tensor);
    const label tensor_offset = (char*)igGrad[0].v_-(char*)(&igGrad[0]);

#if FULL_DEGUG
    Info<<endl<<"!! Init vector RowSubsections for GaussGradf"<<endl<<endl;
#endif
    // calculate the max face number for LDM to include all relative arrays
    static label faceSize4CPE = label
                                (
                                    0.6*64.0*1000.0/
                                    (
										// data on edges
                                        (scalar)vector_size
                                        + (scalar)sizeof(label)
                                        + (scalar)tensor_size
										// data on columns, use edge length
                                        + (scalar)sizeof(scalar)*9.0
                                        + (scalar)tensor_size
                                    )
                                );

    // static RowSubsections object that will not be construct repeatedly
    // when this function is called for more than once
    static SWFoam::RowSubsections secs(fSize, cSize, 64, 4, lPtr, uPtr,
                                      faceSize4CPE, tensor_size);
    // call the C sub-routine
    swGgradfInnerRoutine_vector( Sf_vptr, issf_vptr, igGrad_tptr, lPtr, uPtr,
                                 fSize, cSize, vector_size, vector_offset,
                                 tensor_size, tensor_offset,
                                 secs.getSubsections(), secs.getSecNumInSeg(),
                                 secs.getColRoundNum());
}

template<>
void swGaussGradDivide<vector>
(Field<vector>& igGrad,
 const volScalarField::DimensionedInternalField& volume)
{
    scalar* igGrad_ptr  = (scalar*) (igGrad.data());
    const label cSize  = igGrad.size(); // cell field size

    const scalar* vol_ptr = volume.begin();

    swGgradfDivide_host(igGrad_ptr, vol_ptr, cSize, 3);

}

template<>
void swGaussGradDivide<tensor>
(Field<tensor>& igGrad,
 const volScalarField::DimensionedInternalField& volume)
{
    scalar* igGrad_ptr  = (scalar*) (igGrad.data());
    const label cSize  = igGrad.size(); // cell field size

    const scalar* vol_ptr = volume.begin();

    swGgradfDivide_host(igGrad_ptr, vol_ptr, cSize, 9);
}

} // namespace Foam

