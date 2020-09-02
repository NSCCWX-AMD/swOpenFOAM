#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(interpolation);
define_e2v_slaveFunPtr(interpolation);

#define swInterpolation_host(cOpt, BEDPtr,  FEDPtr, vertexPtr, \
			edgeNum ,vertexNum, dims) \
{ \
	Arrays backEdgeData,frontEdgeData,selfConnData,vertexData; \
	if(!BEDPtr) \
	{ \
		constructEmptyArray(backEdgeData); \
	} else \
	{ \
		constructSingleArray(backEdgeData, dims, (edgeNum) ,COPYIN, \
					BEDPtr); \
	} \
	if(!FEDPtr) \
	{ \
		constructEmptyArray(frontEdgeData); \
	} else \
	{ \
		constructSingleArray(frontEdgeData, dims, (edgeNum) ,COPYOUT, \
					FEDPtr); \
	} \
	constructEmptyArray(selfConnData); \
	if(!vertexPtr) \
	{ \
		constructEmptyArray(vertexData); \
	} else \
	{ \
		constructSingleArray(vertexData, dims, (vertexNum),COPYIN, \
						vertexPtr); \
	} \
	FieldData data \
		= {&backEdgeData, &frontEdgeData, &selfConnData, &vertexData}; \
	(cOpt)->fun_slave = slave_interpolation; \
	(cOpt)->fun_host  = interpolation; \
	(cOpt)->data      = &data; \
}

#ifdef __cplusplus
}
#endif

#endif
