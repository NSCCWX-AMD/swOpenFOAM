#ifndef CALCULATEUVWFLUX_H
#define CALCULATEUVWFLUX_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swCalculateUvwFlux);
define_e2v_slaveFunPtr(swCalculateUvwFlux);

#define swCalculateUvwFlux_data(massFlux, u, v, w, Su, Sv, Sw, \
			edgeNum, vertexNum) \
{ \
	int i,j; \
	u  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	v  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	w  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	Su = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	Sv = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	Sw = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* phi[6] = {u, v, w, Su, Sv, Sw} ;\
	for(i=0;i<6;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			phi[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	massFlux = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	for(j=0;j<edgeNum;j++) \
	{ \
		massFlux[j] = (swFloat)(2*j+1)/(2*j+2); \
	} \
} 

#define swCalculateUvwFlux_master(cOpt, massFlux, u, v, w, Su, Sv, Sw, \
			edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	swFloat* phi[3] = {u, v, w}; \
	swFloat* S[3] = {Su, Sv, Sw}; \
	Arrays backEdgeData[3], frontEdgeData[3]; \
	Arrays selfConnData[3], vertexData[3]; \
	FieldData data[3]; \
	Arrays backEdgeData_p[3], frontEdgeData_p[3]; \
	Arrays selfConnData_p[3], vertexData_p[3]; \
	FieldData data_p[3]; \
	for(iArray=0;iArray<3;iArray++) \
	{ \
		constructEmptyArray(backEdgeData[iArray]); \
		constructEmptyArray(selfConnData[iArray]); \
		data[iArray].backEdgeData  = &backEdgeData[iArray]; \
		data[iArray].frontEdgeData = &frontEdgeData[iArray]; \
		data[iArray].selfConnData  = &selfConnData[iArray]; \
		data[iArray].vertexData    = &vertexData[iArray]; \
		cOpt[iArray].data = &data[iArray]; \
		\
		constructEmptyArray(backEdgeData_p[iArray]); \
		constructEmptyArray(selfConnData_p[iArray]); \
		data_p[iArray].backEdgeData  = &backEdgeData_p[iArray]; \
		data_p[iArray].frontEdgeData = &frontEdgeData_p[iArray]; \
		data_p[iArray].selfConnData  = &selfConnData_p[iArray]; \
		data_p[iArray].vertexData    = &vertexData_p[iArray]; \
		cOpt[iArray].data_p = &data_p[iArray]; \
	} \
	/* operator 0 */ \
	opt = 0; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, u); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				massFlux); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_host  = swCalculateUvwFlux; \
	cOpt[opt].fun_slave = slave_swCalculateUvwFlux; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, v); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_host  = swCalculateUvwFlux; \
	cOpt[opt].fun_slave = slave_swCalculateUvwFlux; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, w); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_host  = swCalculateUvwFlux; \
	cOpt[opt].fun_slave = slave_swCalculateUvwFlux; \
	\
}

#define swCalculateUvwFlux_test() \
{ \
    int optNum = 3; \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *u, *v, *w, *Su, *Sv, *Sw, *massFlux; \
	swFloat *u_m, *v_m, *w_m, *Su_m, *Sv_m, *Sw_m, *massFlux_m; \
	swFloat gamblend = 1.1234566; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	swFloat paras[1] = {gamblend}; \
	constructSingleArray(paraData, 1, 1, COPYIN, &paras[0]); \
	swCalculateUvwFlux_data(massFlux, u, v, w, Su, Sv, Sw, \
				edgeNum, vertexNum); \
    swCalculateUvwFlux_data(massFlux_m, u_m, v_m, w_m, Su_m, Sv_m, Sw_m, \
				edgeNum, vertexNum); \
	swCalculateUvwFlux_master(cOpt_m, massFlux_m, u_m, v_m, w_m, \
				Su_m, Sv_m, Sw_m, edgeNum, vertexNum); \
	\
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
	swCalculateUvwFlux_host(&mlbIter, massFlux, u, v, w, Su, Sv, Sw, \
				gamblend, edgeNum, vertexNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(Su_m, Su, vertexNum*dims); \
	checkResult(Sv_m, Sv, vertexNum*dims); \
	checkResult(Sw_m, Sw, vertexNum*dims); \
}

#ifdef __cplusplus
}
#endif

#endif
