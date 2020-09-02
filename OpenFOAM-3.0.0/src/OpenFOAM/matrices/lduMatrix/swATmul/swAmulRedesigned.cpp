#include "GAMGAgglomeration.H"
#include "lduMatrix.H"
#include "labelField.H"
#include "primitiveFields.H"
#include "LUscalarMatrix.H"
#include "GAMGSolver.H"
#include "swGAMGSolver.hpp"
#include <mpi.h>
#include <vector>

#include "processorFvPatchField.H"
#include "processorGAMGInterfaceField.H"

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

#define INTERFACES_LOOKUP \
scalar *sendBuffer, *recvBuffer;\
label numInterfaces = 0;\
std::vector<label> interfacesUse;\
\
forAll(interfaces, interfaceI)\
{\
	if (interfaces.set(interfaceI))\
	{\
		numInterfaces++;\
		interfacesUse.push_back(interfaceI);\
	}\
}\
\
label *locPosition = new label[numInterfaces+1];\
label *destRank = new label[numInterfaces];\
MPI_Request requestAmul[2*numInterfaces];\
\
locPosition[0] = 0;


#define FVPATCH \
const processorFvPatchField<scalar> *dp = static_cast<const processorFvPatchField<scalar> *>(bp);

#define GAMGPATCH \
const processorGAMGInterfaceField *dp = static_cast<const processorGAMGInterfaceField *>(bp);

#define FVPATCH_FACELLS \
const labelUList& faceCells = dp->patch().faceCells();

#define GAMG_FACELLS \
const labelUList& faceCells = dp->faceCells();

#define INTERFACES_INIT(PATCHPOINTER,FACECELLS) \
for(int i=0; i<numInterfaces; ++i)\
{\
	label interfaceI = interfacesUse[i];\
	const lduInterfaceField *bp = &(interfaces[interfaceI]);\
	PATCHPOINTER\
\
	destRank[i] = dp->neighbProcNo();\
	locPosition[i+1] = dp->size() + locPosition[i];\
}\
\
label bufferSize = locPosition[numInterfaces];\
\
sendBuffer = new scalar[bufferSize];\
recvBuffer = new scalar[bufferSize];\
for(int i=0; i<numInterfaces; ++i)\
{\
	label interfaceI = interfacesUse[i];\
	const lduInterfaceField *bp = &(interfaces[interfaceI]);\
	PATCHPOINTER\
	FACECELLS\
	label locSize = locPosition[i+1] - locPosition[i];\
	for(label facei=0; facei<locSize; facei++)\
    {\
        sendBuffer[facei+locPosition[i]] = psi[faceCells[facei]];\
    }\
}

#define INTERFACES_UPDATE(PATCHPOINTER,FACECELLS) \
for(int i=0; i<numInterfaces; ++i)\
{\
	label interfaceI = interfacesUse[i];\
	const lduInterfaceField *bp = &(interfaces[interfaceI]);\
	PATCHPOINTER\
	FACECELLS\
\
	label locSize = locPosition[i+1] - locPosition[i];\
	for(label facei=0; facei<locSize; facei++)\
    {\
        ApsiPtr[faceCells[facei]] -=\
        	interfaceBouCoeffs[interfaceI][facei] * recvBuffer[facei+locPosition[i]];\
    }\
}

#define ISEND_IRECV \
for(int i=0; i<numInterfaces; ++i)\
{\
	label locSize = locPosition[i+1] - locPosition[i];\
	MPI_Isend\
	(\
		&(sendBuffer[locPosition[i]]),\
		locSize,\
		MPI_DOUBLE,\
		destRank[i],\
		1,\
		MPI_COMM_WORLD,\
		&requestAmul[i]\
	);\
	MPI_Irecv\
	(\
		&(recvBuffer[locPosition[i]]),\
		locSize,\
		MPI_DOUBLE,\
		destRank[i],\
		1,\
		MPI_COMM_WORLD,\
		&requestAmul[i + numInterfaces]\
	);\
}


#define SWCOMPUTE \
amul_parameter.psiPtr = (SCALAR*)(&psiPtr[0]);\
amul_parameter.ApsiPtr =(SCALAR*)(&ApsiPtr[0]);\
for(int i = 0 ; i < CORE_SIZE ; i++) {\
    matrix_translate[i].psiPtr = amul_parameter.psiPtr;\
    matrix_translate[i].ApsiPtr = amul_parameter.ApsiPtr;\
}\
amul_host(&amul_parameter, matrix_translate);


#define HOST_COMPUTE \
const scalar* const __restrict__ diagPtr = matrix.diag().begin();\
\
const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();\
const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();\
\
const scalar* const __restrict__ upperPtr = matrix.upper().begin();\
const scalar* const __restrict__ lowerPtr = matrix.lower().begin();\
\
register const label nCells = matrix.diag().size();\
\
for (register label cell=0; cell<nCells; cell++)\
{\
    ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];\
}\
\
register const label nFaces = matrix.upper().size();\
\
for (register label face=0; face<nFaces; face++)\
{\
    ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];\
    ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];\
}


namespace Foam
{

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
)
{
    scalar* __restrict__ ApsiPtr = Apsi.begin();
    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();

    if(UPstream::parRun())
    {
        INTERFACES_LOOKUP
        INTERFACES_INIT(FVPATCH,FVPATCH_FACELLS)
        ISEND_IRECV

        SWCOMPUTE
        // HOST_COMPUTE

        MPI_Waitall(2*numInterfaces, &requestAmul[0], MPI_STATUSES_IGNORE);

        INTERFACES_UPDATE(FVPATCH,FVPATCH_FACELLS)

        tpsi.clear();
        delete [] sendBuffer;
        delete [] recvBuffer;
        delete [] locPosition;
        delete [] destRank;
    }
    else
    {
        SWCOMPUTE
        // HOST_COMPUTE
    }
}


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
)
{
    scalar* __restrict__ ApsiPtr = Apsi.begin();
    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();

    if(UPstream::parRun())
    {
        INTERFACES_LOOKUP
        INTERFACES_INIT(GAMGPATCH,GAMG_FACELLS)
        ISEND_IRECV

        SWCOMPUTE
        // HOST_COMPUTE

        MPI_Waitall(2*numInterfaces, &requestAmul[0], MPI_STATUSES_IGNORE);

        INTERFACES_UPDATE(GAMGPATCH,GAMG_FACELLS)

        tpsi.clear();
        delete [] sendBuffer;
        delete [] recvBuffer;
        delete [] locPosition;
        delete [] destRank;
    }
    else
    {
        SWCOMPUTE
        // HOST_COMPUTE
    }
}


void GAMGAmulRedesigned
(
    const lduMatrix & matrix,
    scalarField& Apsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
)
{
    scalar* __restrict__ ApsiPtr = Apsi.begin();
    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();

    if(UPstream::parRun())
    {
        INTERFACES_LOOKUP
        INTERFACES_INIT(GAMGPATCH,GAMG_FACELLS)
        ISEND_IRECV
        HOST_COMPUTE
        MPI_Waitall(2*numInterfaces, &requestAmul[0], MPI_STATUSES_IGNORE);
        INTERFACES_UPDATE(GAMGPATCH,GAMG_FACELLS)
        tpsi.clear();
        delete [] sendBuffer;
        delete [] recvBuffer;
        delete [] locPosition;
        delete [] destRank;
    }
    else
    {
        HOST_COMPUTE
    }
}


void initAndUpdate
(
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const tmp<scalarField>& tpsi,
    scalarField& Apsi,
    const direction cmpt
)
{
    scalar* __restrict__ ApsiPtr = Apsi.begin();
    const scalarField& psi = tpsi();

    if(UPstream::parRun())
    {
        INTERFACES_LOOKUP
        INTERFACES_INIT(FVPATCH,FVPATCH_FACELLS)
        ISEND_IRECV

        MPI_Waitall(2*numInterfaces, &requestAmul[0], MPI_STATUSES_IGNORE);
        INTERFACES_UPDATE(FVPATCH,FVPATCH_FACELLS)

        tpsi.clear();
        delete [] sendBuffer;
        delete [] recvBuffer;
        delete [] locPosition;
        delete [] destRank;
    }
}

} // end namespace Foam

