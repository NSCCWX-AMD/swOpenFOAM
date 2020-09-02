#ifndef SEPARATE_INTEGRATE_H
#define SEPARATE_INTEGRATE_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swSeparateIntegrate_su);
define_e2v_slaveFunPtr(swSeparateIntegrate_su);

#define swSeparateIntegrate_data(fu, su, edgeNum, vertexNum) \
{ \
	int i,j; \
	fu        = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* edge[1] = {fu}; \
	for(i=0;i<1;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	su      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* cell[1] = {su}; \
	for(i=0;i<1;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			cell[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
} 

#define swSeparateIntegrate_master(cOpt, fu, su, edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	Arrays backEdgeData[1], frontEdgeData[1]; \
	Arrays selfConnData[1], vertexData[1]; \
	FieldData data[1]; \
	Arrays backEdgeData_p[1], frontEdgeData_p[1]; \
	Arrays selfConnData_p[1], vertexData_p[1]; \
	FieldData data_p[1]; \
	for(iArray=0;iArray<1;iArray++) \
	{ \
		data[iArray].backEdgeData    = &backEdgeData[iArray]; \
		data[iArray].frontEdgeData   = &frontEdgeData[iArray]; \
		data[iArray].selfConnData    = &selfConnData[iArray]; \
		data[iArray].vertexData      = &vertexData[iArray]; \
		cOpt[iArray].data = &data[iArray]; \
		data_p[iArray].backEdgeData  = &backEdgeData_p[iArray]; \
		data_p[iArray].frontEdgeData = &frontEdgeData_p[iArray]; \
		data_p[iArray].selfConnData  = &selfConnData_p[iArray]; \
		data_p[iArray].vertexData    = &vertexData_p[iArray]; \
		cOpt[iArray].data_p = &data_p[iArray]; \
	} \
	/* operator 0 */ \
	opt = 0; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, fu); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, su); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swSeparateIntegrate_su; \
	cOpt[opt].fun_host  = swSeparateIntegrate_su; \
	\
}

#define swSeparateIntegrate_test() \
{ \
    int optNum = 1; \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *su, *fu; \
	swFloat *su_m, *fu_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	swSeparateIntegrate_data(fu, su, edgeNum, vertexNum); \
	swSeparateIntegrate_data(fu_m, su_m, edgeNum, vertexNum); \
	swSeparateIntegrate_master(cOpt_m, fu_m, su_m, edgeNum, vertexNum); \
	printf("test\n"); \
	getTime(time1); \
	for(int iOpt=0;iOpt<optNum;iOpt++) \
	{ \
		backEdgeData_m  = *(cOpt_m[iOpt].data_p->backEdgeData); \
		frontEdgeData_m = *(cOpt_m[iOpt].data_p->frontEdgeData); \
		selfConnData_m  = *(cOpt_m[iOpt].data_p->selfConnData); \
		vertexData_m    = *(cOpt_m[iOpt].data_p->vertexData); \
		cOpt_m[iOpt].fun_host(&backEdgeData_m, &frontEdgeData_m, \
					&selfConnData_m, &vertexData_m, &paraData, \
					rowAddr,colAddr, cOpt_m[iOpt].data); \
	} \
    getTime(time2); \
	printf("CPU Processor Time: %f us\n", (time2-time1)*1000000); \
	CG_init(); \
    getTime(time1); \
	swSeparateIntegrate_host(&mlbIter, fu, su, edgeNum, vertexNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(su_m, su, vertexNum); \
}

#ifdef __cplusplus
}
#endif

#endif
