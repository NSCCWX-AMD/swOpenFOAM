#ifndef ITERATOR_H
#define ITERATOR_H

#include "swMacro.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C"{
#endif

//macros to define core function pointer
//	use as:
//		define_e2v_hostFunPtr( spMV )
//		{
//			Here do the core computing;
//			Parameters is similiar as edge2VertexIteration.
//		}
#define define_e2v_FunPtr( funname ) \
void funname (Arrays* backEdgeData, Arrays* frontEdgeData, \
			Arrays* selfConnData, Arrays* vertexData, Arrays* paraData, \
			swInt* startVertices, swInt* endVertices, FieldData *data)

#define define_e2v_slaveFunPtr( funname ) \
void slave_##funname (Arrays* backEdgeData, Arrays* frontEdgeData, \
			Arrays* selfConnData, Arrays* vertexData, Arrays *paraData, \
			swInt* startVertices, swInt* endVertices, FieldData *data)

#define define_array_FunPtr( funname ) \
void funname (Arrays* backEdgeData, Arrays* frontEdgeData, \
			Arrays* selfConnData, Arrays* vertexData, Arrays* paraData, \
			FieldData *data)

#define define_array_slaveFunPtr( funname ) \
void slave_##funname (Arrays* backEdgeData, Arrays* frontEdgeData, \
			Arrays* selfConnData, Arrays* vertexData, Arrays *paraData, \
			FieldData *data)

#define define_v2e_hostFunPtr( funname ) \
void funname (Arrays* neighbourData, Arrays* vertexData, \
			swInt* accuEdgeNumbers, swInt* neighbourVetices)

#define define_v2e_slaveFunPtr( funname ) \
void funname (Arrays* neighbourData, Arrays* vertexData, \
			swInt* accuEdgeNumbers, swInt* neighbourVetices)

// the indicator to determine data copy action
#define COPYIN 0 
#define COPYOUT 1 
#define COPYINOUT 2
#define UPDATED 3


// parameter struct for interation
typedef struct
{
	// float data array holder 
	swFloat** floatArrays;
	// integer data array holder 
	swInt** intArrays;
	// float data array dimensions 
	swInt*    fArrayDims;
	// integer data array dimensions 
	swInt*    iArrayDims;
	// float data array action 
	swInt*    fArrayInOut;
	// integer data array action
	swInt*    iArrayInOut;
	// float data array number 
	swInt    fArrayNum;
	// integer data array number 
	swInt    iArrayNum;
	// float data array sizes 
	swInt    fArraySizes;
	// integer data array sizes
	swInt    iArraySizes;
} Arrays;

typedef struct
{
	Arrays *backEdgeData;
	Arrays *frontEdgeData;
	Arrays *selfConnData;
	Arrays *vertexData;
} FieldData;

//function pointer
typedef void (* e2v_hostFunPtr)(Arrays* backEdgeData, Arrays* frontEdgeData,
			Arrays* selfConnData, Arrays* vertexData, Arrays* paraData, 
			swInt* startVertices, swInt* endVertices, FieldData *data);
typedef void (* e2v_slaveFunPtr)(Arrays* backEdgeData,Arrays* frontEdgeData,
			Arrays* selfConnData, Arrays* vertexData, Arrays* paraData, 
			swInt* startVertices, swInt* endVertices, FieldData *data);
typedef void (* array_hostFunPtr)(Arrays* backEdgeData, 
			Arrays* frontEdgeData, Arrays* selfConnData,Arrays* vertexData,
		   	Arrays* paraData, FieldData *data);
typedef void (* array_slaveFunPtr)(Arrays* backEdgeData,
			Arrays* frontEdgeData, Arrays* selfConnData,Arrays* vertexData,
		   	Arrays* paraData, FieldData *data);
typedef void (* v2e_hostFunPtr) (Arrays* neighbourData, Arrays* vertexData,
			swInt* accuEdgeNumbers, swInt* neighbourVetices);
typedef void (* v2e_slaveFunPtr) (Arrays* neighbourData, Arrays* vertexData,
			swInt* accuEdgeNumbers, swInt* neighbourVetices);

// function pointer binded with arrays
typedef struct
{
	e2v_hostFunPtr fun_host;
	e2v_slaveFunPtr fun_slave;
	FieldData *data;
	FieldData *data_p;
} coupledOperator;

typedef struct
{
	array_hostFunPtr fun_host;
	array_slaveFunPtr fun_slave;
	FieldData *data;
} coupledArrayOperator;


/*************** Macro to construct empty arrays ***************/
#define constructEmptyArray( arrays ) \
{ \
	(arrays).iArraySizes = 0; \
	(arrays).iArrayNum = 0; \
	(arrays).iArrayInOut = NULL; \
	(arrays).iArrayDims = NULL; \
	(arrays).intArrays = NULL;\
	\
	\
	(arrays).fArraySizes = 0; \
	(arrays).fArrayNum = 0; \
	(arrays).fArrayInOut = NULL; \
	(arrays).fArrayDims = NULL; \
	(arrays).floatArrays = NULL;\
}

// Construct Arrays from source arrays but different size
#define constructFromArrays(srcArrays,dstArrays, size) \
{\
	int i; \
	(dstArrays)->fArraySizes = (srcArrays)->fArraySizes; \
	(dstArrays)->fArrayNum   = (srcArrays)->fArrayNum; \
	(dstArrays)->fArrayInOut = NEW(swInt,    (dstArrays->fArrayNum)); \
	(dstArrays)->fArrayDims  = NEW(swInt,    (dstArrays->fArrayNum)); \
	(dstArrays)->floatArrays = NEW(swFloat*, (dstArrays->fArrayNum)); \
	for(i=0;i<(dstArrays)->fArrayNum;i++) \
	{ \
		(dstArrays)->fArrayInOut[i] = (srcArrays)->fArrayInOut[i]; \
		(dstArrays)->fArrayDims[i]  = (srcArrays)->fArrayDims[i]; \
		(dstArrays)->floatArrays[i] \
		= NEW(swFloat,(dstArrays)->fArraySizes*dstArrays->fArrayDims[i]); \
	} \
	\
	\
	(dstArrays)->iArraySizes = (srcArrays)->iArraySizes; \
	(dstArrays)->iArrayNum   = (srcArrays)->iArrayNum; \
	(dstArrays)->iArrayInOut = NEW(swInt, (dstArrays->iArrayNum)); \
	(dstArrays)->iArrayDims  = NEW(swInt, (dstArrays->iArrayNum)); \
	(dstArrays)->intArrays   = NEW(swInt*, (dstArrays->iArrayNum)); \
	for(i=0;i<(dstArrays)->iArrayNum;i++) \
	{ \
		(dstArrays)->iArrayInOut[i] = (srcArrays)->iArrayInOut[i]; \
		(dstArrays)->iArrayDims[i]  = (srcArrays)->iArrayDims[i]; \
		(dstArrays)->intArrays[i] \
		= NEW(swInt,(dstArrays)->iArraySizes*dstArrays->iArrayDims[i]); \
	} \
}

// release the memory space using delete
#define deleteArrays(Arrays) \
{ \
	int i; \
	for(i=0;i<(Arrays)->fArrayNum;i++) \
	{ \
		delete((Arrays)->floatArrays[i]); \
	} \
	delete((Arrays)->fArrayInOut); \
	delete((Arrays)->fArrayDims); \
	delete((Arrays)->floatArrays); \
	\
	for(i=0;i<(Arrays)->iArrayNum;i++) \
	{ \
		delete((Arrays)->intArrays[i]); \
	} \
	delete((Arrays)->iArrayInOut); \
	delete((Arrays)->iArrayDims); \
	delete((Arrays)->intArrays); \
}

/****************** Macros handle float arrays ******************/
// address copy
#define constructSingleArray( arrays, dim, size, io, pointer) \
{ \
	(arrays).iArraySizes = 0; \
	(arrays).iArrayNum = 0; \
	(arrays).iArrayInOut = NULL; \
	(arrays).iArrayDims = NULL; \
	(arrays).intArrays = NULL;\
	\
	\
	(arrays).fArraySizes = size; \
	(arrays).fArrayNum = 1; \
	\
	(arrays).fArrayInOut = NEW(swInt, 1); \
	*((arrays).fArrayInOut) = io; \
	\
	(arrays).fArrayDims = NEW(swInt, 1); \
	*((arrays).fArrayDims) = dim; \
	\
	(arrays).floatArrays = NEW(swFloat*, 1);\
	*(arrays).floatArrays = pointer; \
}

// address copy: arrays to arrays
#define constructSingleArrays( arrays, srcArrays) \
{ \
	(arrays)->iArraySizes = 0; \
	(arrays)->iArrayNum = 0; \
	(arrays)->iArrayInOut = NULL; \
	(arrays)->iArrayDims = NULL; \
	(arrays)->intArrays = NULL;\
	\
	\
	(arrays)->fArraySizes = (srcArrays)->fArraySizes; \
	(arrays)->fArrayNum = (srcArrays)->fArrayNum; \
	\
	(arrays)->fArrayInOut = (srcArrays)->fArrayInOut; \
	\
	(arrays)->fArrayDims = (srcArrays)->fArrayDims; \
	\
	(arrays)->floatArrays = NEW(swFloat*, (arrays)->fArrayNum);\
	int iArray; \
	for(iArray=0;iArray<(arrays)->fArrayNum;iArray++) \
	{ \
		(arrays)->floatArrays[iArray] = (srcArrays)->floatArrays[iArray]; \
	} \
}

// address copy
#define addSingleArray( arrays, dim, size, io, pointer) \
{ \
	if( (arrays).fArraySizes != size) \
	{ \
		dumpError("can not add array with different length!\n"); \
		exit(-1); \
	} \
	\
	(arrays).fArrayNum++; \
	swInt fArrayNum = (arrays).fArrayNum; \
	\
	RESIZE(swInt, (arrays).fArrayInOut, fArrayNum-1, fArrayNum); \
	((arrays).fArrayInOut)[fArrayNum-1] = io; \
	\
	RESIZE(swInt, (arrays).fArrayDims, fArrayNum-1, fArrayNum); \
	((arrays).fArrayDims)[fArrayNum-1] = dim; \
	\
	RESIZE(swFloat*, (arrays).floatArrays, fArrayNum-1, fArrayNum);\
	(arrays).floatArrays[fArrayNum-1] = pointer;\
}

// deep copy
#define copySingleArray( arrays, dim, size, io, pointer) \
{ \
	if( (arrays).fArraySizes != (size) ); \
	{ \
		dumpError("can not add array with different length!\n"); \
		exit(-1); \
	} \
	\
	(arrays).fArrayNum++; \
	\
	RESIZE(swInt, (arrays).fArrayInOut, fArrayNum-1, fArrayNum); \
	((arrays).fArrayInOut)[fArrayNum-1] = io; \
	\
	RESIZE(swInt, (arrays).fArrayDims, fArrayNum-1, fArrayNum); \
	((arrays).fArrayDims)[fArrayNum-1] = dim; \
	\
	RESIZE(swFloat*, (arrays).floatArrays, fArrayNum-1, fArrayNum);\
	(arrays).floatArrays[fArrayNum-1] = NEW( swFloat; size); \
	memcpy( (arrays).floatArrays[fArrayNum-1], pointer, \
				sizeof(swFloat)*(size) ); \
}

// add blank array
#define addBlankArray( arrays, dim, size, io) \
{ \
	if( (arrays).fArraySizes != (size) ); \
	{ \
		dumpError("can not add array with different length!\n"); \
		exit(-1); \
	} \
	\
	(arrays).fArrayNum++; \
	\
	RESIZE(swInt, (arrays).fArrayInOut, fArrayNum-1, fArrayNum); \
	((arrays).fArrayInOut)[fArrayNum-1] = io; \
	\
	RESIZE(swInt, (arrays).fArrayDims, fArrayNum-1, fArrayNum); \
	((arrays).fArrayDims)[fArrayNum-1] = dim; \
	\
	RESIZE(swFloat*, (arrays).floatArrays, fArrayNum-1, fArrayNum);\
	(arrays).floatArrays[fArrayNum-1] = NEW( swFloat; size); \
}


/****************** Macros handle int arrays ******************/
// address copy
#define constructSingleIArray( arrays, dim, size, io, pointer) \
{ \
	(arrays).fArraySizes = 0; \
	(arrays).fArrayNum = 0; \
	(arrays).fArrayInOut = NULL; \
	(arrays).fArrayDims = NULL; \
	(arrays).floatArrays = NULL;\
	\
	\
	(arrays).iArraySizes = size; \
	(arrays).iArrayNum = 1; \
	\
	(arrays).iArrayInOut = NEW(swInt, 1); \
	*((arrays).iArrayInOut) = io; \
	\
	(arrays).iArrayDims = NEW(swInt, 1); \
	*((arrays).iArrayDims) = dim; \
	\
	(arrays).intArrays = NEW(swInt*, 1);\
	*(arrays).intArrays = pointer;\
}

// address copy
#define addSingleIArray( arrays, dim, size, io, pointer) \
{ \
	if( (arrays).iArraySizes != size); \
	{ \
		dumpError("can not add array with different length!\n"); \
		exit(-1); \
	} \
	\
	(arrays).iArrayNum++; \
	\
	RESIZE(swInt, (arrays).iArrayInOut, iArrayNum-1, iArrayNum); \
	((arrays).iArrayInOut)[iArrayNum-1] = io; \
	\
	RESIZE(swInt, (arrays).iArrayDims, iArrayNum-1, iArrayNum); \
	((arrays).iArrayDims)[iArrayNum-1] = dim; \
	\
	RESIZE(swInt*, (arrays).intArrays, iArrayNum-1, iArrayNum);\
	(arrays).intArrays[iArrayNum-1] = pointer;\
}

// deep copy
#define copySingleIArray( arrays, dim, size, io, pointer) \
{ \
	if( (arrays).iArraySizes != (size) ); \
	{ \
		dumpError("can not add array with different length!\n"); \
		exit(-1); \
	} \
	\
	(arrays).iArrayNum++; \
	\
	RESIZE(swInt, (arrays).iArrayInOut, iArrayNum-1, iArrayNum); \
	((arrays).iArrayInOut)[iArrayNum-1] = io; \
	\
	RESIZE(swInt, (arrays).iArrayDims, iArrayNum-1, iArrayNum); \
	((arrays).iArrayDims)[iArrayNum-1] = dim; \
	\
	RESIZE(swFloat*, (arrays).intArrays, iArrayNum-1, iArrayNum);\
	(arrays).intArrays[iArrayNum-1] = NEW( swInt; size); \
	memcpy( (arrays).iArrays[iArrayNum-1], pointer, \
				sizeof(swInt)*(size) ); \
}

// add blank array
#define addBlankIArray( arrays, dim, size, io) \
{ \
	if( (arrays).iArraySizes != (size) ); \
	{ \
		dumpError("can not add array with different length!\n"); \
		exit(-1); \
	} \
	\
	(arrays).iArrayNum++; \
	\
	RESIZE(swInt, (arrays).iArrayInOut, iArrayNum-1, iArrayNum); \
	((arrays).iArrayInOut)[iArrayNum-1] = io; \
	\
	RESIZE(swInt, (arrays).iArrayDims, iArrayNum-1, iArrayNum); \
	((arrays).iArrayDims)[iArrayNum-1] = dim; \
	\
	RESIZE(swFloat*, (arrays).intArrays, iArrayNum-1, iArrayNum);\
	(arrays).intArrays[iArrayNum-1] = NEW( swInt; size); \
}


/****************** Macros to destroy arrays ******************/
// shalow free
#define destroyArray( arrays ) \
{ \
	DELETE( (arrays).fArrayInOut ); \
	DELETE( (arrays).fArrayDims  ); \
	DELETE( (arrays).floatArrays ); \
	DELETE( (arrays).iArrayInOut ); \
	DELETE( (arrays).iArrayDims  ); \
	DELETE( (arrays).intArrays ); \
}

// deep memory free
#define deleteSingleArray( arrays ) \
{ \
	size_t i = (array)->fArrayNum;\
	while(i--){ DELETE((array)->floatArrays[i])}; \
	DELETE(floatArrays); \
	DELETE(fArrayDims); \
	DELETE(fArrayInOut); \
	\
	i = (array)->iArrayNum; \
	while(i--){ DELETE((array)->intArrays[i])}; \
	DELETE(intArrays); \
	DELETE(iArrayDims); \
	DELETE(iArrayInOut); \
}

// array data accessor
#define accessArray( array, id ) \
( (array)->floatArrays[id] )

#define getArraySize( array ) \
( (array)->fArraySizes )

#define getArrayNumber( array ) \
( (array)->fArrayNum )

#define getArrayDims( array, id) \
( (array)->fArrayDims[id] )

#define getArrayInOut( array, id ) \
( (array)->fArrayInOut[id] )

#define accessIArray( array, id ) \
( (array)->intArrays[id] )

#define getIArraySize( array ) \
( (array)->iArraySizes )

#define getIArrayNumber( array ) \
( (array)->iArrayNum )

#define getIArrayDims( array ) \
( (array)->iArrayDims )

#define getIArrayInOut( array ) \
( (array)->iArrayInOut )

#ifdef __cplusplus
}
#endif

#endif // ITERATOR_H
