#ifndef SPMV_H
#define SPMV_H
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(swSpMV);
define_e2v_slaveFunPtr(swSpMV);

#define swSpMV_data(lowerPtr, upperPtr, bPtr, xPtr, diagPtr, \
			edgeNum, vertexNum) \
{ \
	int i,j,idim; \
	lowerPtr = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat)); \
	upperPtr = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat)); \
	for(i=0;i<edgeNum;i++) \
	{ \
		for(idim=0;idim<dims;idim++) \
		{ \
			lowerPtr[i*dims+idim] = (swFloat)(i+1)/(i+2); \
			upperPtr[i*dims+idim] = (swFloat)(i+2)/(i+1); \
		} \
	} \
	xPtr = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat)); \
	bPtr = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat)); \
	diagPtr = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat)); \
	for(i=0;i<vertexNum;i++) \
	{ \
		for(idim=0;idim<dims;idim++) \
		{ \
			bPtr[i*dims+idim] = (swFloat)(2*i+1)/(2*i+2); \
			xPtr[i*dims+idim] = (swFloat)(2*i+2)/(2*i+3); \
			diagPtr[i*dims+idim] = (swFloat)(i+2)/(i+3); \
		} \
	} \
}

#define swSpMV_v2e_data(edge, bPtr, xPtr, diagPtr, \
			edgeNum, vertexNum) \
{ \
	int i,j,idim; \
	edge = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat)); \
	for(i=0;i<edgeNum;i++) \
	{ \
		for(idim=0;idim<dims;idim++) \
		{ \
			edge[i*dims+idim] = (swFloat)(i+1)/(i+2); \
		} \
	} \
	xPtr = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat)); \
	bPtr = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat)); \
	diagPtr = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat)); \
	for(i=0;i<vertexNum;i++) \
	{ \
		for(idim=0;idim<dims;idim++) \
		{ \
			bPtr[i*dims+idim] = (swFloat)(2*i+1)/(2*i+2); \
			xPtr[i*dims+idim] = (swFloat)(2*i+2)/(2*i+3); \
			diagPtr[i*dims+idim] = (swFloat)(i+2)/(i+3); \
		} \
	} \
}
#define swSpMV_master(cOpt, lowerPtr, upperPtr, bPtr, xPtr, diagPtr, \
			edgeNum ,vertexNum) \
{ \
	Arrays backEdgeData,frontEdgeData,selfConnData,vertexData; \
	constructSingleArray(backEdgeData, dims, (edgeNum) ,COPYIN, \
				lowerPtr); \
	constructSingleArray(frontEdgeData, dims, (edgeNum) ,COPYIN, \
				upperPtr); \
	constructSingleArray(selfConnData, dims, (vertexNum) ,COPYIN, \
				diagPtr); \
	constructSingleArray(vertexData, dims, (vertexNum),COPYIN,xPtr); \
	addSingleArray(vertexData, dims, (vertexNum), COPYOUT, bPtr); \
	FieldData data \
		= {&backEdgeData, &frontEdgeData, &selfConnData, &vertexData}; \
	cOpt[0].fun_slave = slave_swSpMV; \
	cOpt[0].fun_host  = swSpMV; \
	cOpt[0].data      = &data; \
}

#define swSpMV_v2e_master(cOpt, edge, bPtr, xPtr, diagPtr, \
			edgeNum ,vertexNum) \
{ \
	Arrays backEdgeData,frontEdgeData,selfConnData,vertexData; \
	constructEmptyArray(backEdgeData); \
	constructSingleArray(frontEdgeData, dims, (edgeNum) ,COPYIN, edge); \
	constructSingleArray(selfConnData, dims, (vertexNum) ,COPYIN, \
				diagPtr); \
	constructSingleArray(vertexData, dims, (vertexNum),COPYIN,xPtr); \
	addSingleArray(vertexData, dims, (vertexNum), COPYOUT, bPtr); \
	FieldData data \
		= {&backEdgeData, &frontEdgeData, &selfConnData, &vertexData}; \
	cOpt[0].fun_slave = slave_swSpMV; \
	cOpt[0].fun_host  = swSpMV; \
	cOpt[0].data      = &data; \
}
#define swSpMV_test(pattern) \
{ \
	int optNum =1; \
	coupledOperator *cOpt = (coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	coupledOperator *cOpt_m = (coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	if(pattern==E2V) \
	{ \
		swFloat *lowerPtr, *upperPtr, *xPtr, *bPtr, *diagPtr; \
		swFloat *lowerPtr_m, *upperPtr_m, *xPtr_m, *bPtr_m, *diagPtr_m; \
		Arrays paraData; \
		constructEmptyArray(paraData); \
		swSpMV_data(lowerPtr, upperPtr, bPtr, xPtr, diagPtr, \
					edgeNum, vertexNum); \
		swSpMV_data(lowerPtr_m, upperPtr_m, bPtr_m, xPtr_m, diagPtr_m, \
					edgeNum, vertexNum); \
		swSpMV_master(cOpt_m, lowerPtr_m, upperPtr_m, bPtr_m, xPtr_m, \
					diagPtr_m, edgeNum, vertexNum); \
		\
		printf("test\n"); \
	    getTime(time1); \
		for(int iOpt=0;iOpt<optNum;iOpt++) \
		{ \
			cOpt_m[iOpt].fun_host(NULL, NULL, NULL, NULL, NULL, \
						startVertices, endVertices, cOpt_m[iOpt].data); \
		} \
		getTime(time2); \
		printf("CPU Processor Time: %f us\n", (time2-time1)*1000000); \
		CG_init(); \
		getTime(time1); \
		swSpMV_e2v_host(&iter, lowerPtr, upperPtr, bPtr, xPtr, diagPtr, \
					dims, edgeNum, vertexNum); \
		getTime(time2); \
		printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
		CG_halt(); \
		\
		checkResult(bPtr_m, bPtr, vertexNum*dims); \
	} else if(pattern==V2E) \
	{ \
		swFloat *edge, *xPtr, *bPtr, *diagPtr; \
		swFloat *edge_m, *xPtr_m, *bPtr_m, *diagPtr_m; \
		Arrays paraData; \
		constructEmptyArray(paraData); \
		swSpMV_v2e_data(edge, bPtr, xPtr, diagPtr, edgeNum, vertexNum); \
		swSpMV_v2e_data(edge_m, bPtr_m, xPtr_m, diagPtr_m, edgeNum, \
					vertexNum); \
		swSpMV_v2e_master(cOpt_m, edge_m, bPtr_m, xPtr_m, diagPtr_m, \
					edgeNum, vertexNum); \
		\
		printf("test\n"); \
	    getTime(time1); \
		for(int iOpt=0;iOpt<optNum;iOpt++) \
		{ \
			cOpt_m[iOpt].fun_host(NULL, NULL, NULL, NULL, NULL, \
						startVertices, endVertices, cOpt_m[iOpt].data); \
		} \
		getTime(time2); \
		printf("CPU Processor Time: %f us\n", (time2-time1)*1000000); \
		CG_init(); \
		getTime(time1); \
		swSpMV_v2e_host(&iter, edge, bPtr, xPtr, diagPtr, \
					dims, edgeNum, vertexNum); \
		getTime(time2); \
		printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
		CG_halt(); \
		\
		checkResult(bPtr_m, bPtr, vertexNum*dims); \
	} \
}

#ifdef __cplusplus
}
#endif

#endif
