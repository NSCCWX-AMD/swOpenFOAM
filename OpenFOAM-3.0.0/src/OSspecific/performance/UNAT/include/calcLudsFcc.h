#ifndef CALCLUDSFCC_H
#define CACLLUDSFCC_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swCalcLudsFcc);
define_e2v_slaveFunPtr(swCalcLudsFcc);
define_e2v_FunPtr(swCalcLudsFcc_rface);
define_e2v_slaveFunPtr(swCalcLudsFcc_rface);
define_e2v_FunPtr(swCalcLudsFcc_rface_CDS);
define_e2v_slaveFunPtr(swCalcLudsFcc_rface_CDS);

#define swCalcLudsFcc_data(massFlux, facex_x, facex_y, facex_z, \
			cellx_x, cellx_y, cellx_z, fccx, fccy, fccz, rface0, rface1, \
			facn, edgeNum, vertexNum) \
{ \
	int i,j; \
	massFlux = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	facex_x  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	facex_y  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	facex_z  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fccx     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fccy     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fccz     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface0   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	facn     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* phi[10] = {massFlux, facex_x, facex_y, facex_z, fccx, fccy, \
		fccz, rface0, rface1, facn}; \
	for(i=0;i<10;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			phi[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	cellx_x = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cellx_y = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cellx_z = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* S[3] = {cellx_x, cellx_y, cellx_z}; \
	for(i=0;i<3;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			S[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
} 

#define swCalcLudsFcc_master(cOpt, massFlux, facex_x, facex_y, facex_z, \
			cellx_x, cellx_y, cellx_z, fccx, fccy, fccz, \
			rface0, rface1, facn, CS, edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	swFloat* facex[3] = {facex_x, facex_y, facex_z}; \
	swFloat* cellx[3] = {cellx_x, cellx_y, cellx_z}; \
	swFloat* fcc[3] = {fccx, fccy, fccz}; \
	Arrays backEdgeData[4], frontEdgeData[4]; \
	Arrays selfConnData[4], vertexData[4]; \
	FieldData data[4]; \
	Arrays backEdgeData_p[4], frontEdgeData_p[4]; \
	Arrays selfConnData_p[4], vertexData_p[4]; \
	FieldData data_p[4]; \
	for(iArray=0;iArray<4;iArray++) \
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
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fccx); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, facex_x); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cellx_x); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				massFlux); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalcLudsFcc; \
	cOpt[opt].fun_host  = swCalcLudsFcc; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fccy); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, facex_y); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cellx_y); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalcLudsFcc; \
	cOpt[opt].fun_host  = swCalcLudsFcc; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fccz); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, facex_z); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cellx_z); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalcLudsFcc; \
	cOpt[opt].fun_host  = swCalcLudsFcc; \
	\
	/* operator 3 */ \
	opt = 3; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface0); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYINOUT, rface1); \
	constructEmptyArray(vertexData[opt]); \
	if(CS==1 || CS==2) \
	{ \
		constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
					massFlux); \
	} else \
	{ \
		constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
					facn); \
	} \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalcLudsFcc_rface; \
	cOpt[opt].fun_host  = swCalcLudsFcc_rface; \
	\
}

#define swCalcLudsFcc_test() \
{ \
    int optNum = 4; \
	coupledOperator *cOpt_m  \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *massFlux, *facex_x, *facex_y, *facex_z; \
	swFloat *cellx_x, *cellx_y, *cellx_z, *fccx, *fccy, *fccz; \
	swFloat *rface0, *rface1, *facn; \
	swFloat *massFlux_m, *facex_x_m, *facex_y_m, *facex_z_m; \
	swFloat *cellx_x_m, *cellx_y_m, *cellx_z_m, *fccx_m, *fccy_m, *fccz_m; \
	swFloat *rface0_m, *rface1_m, *facn_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
   	swFloat gamblend = 1.234567; \
	swFloat CS = 1; \
	swFloat paras[2] = {gamblend, CS}; \
	constructSingleArray(paraData, 1, 2, COPYIN, &paras[0]); \
	swCalcLudsFcc_data(massFlux, facex_x, facex_y, facex_z, cellx_x,  \
				cellx_y, cellx_z, fccx, fccy, fccz, rface0, rface1, facn, \
				edgeNum, vertexNum); \
	swCalcLudsFcc_data(massFlux_m, facex_x_m, facex_y_m, facex_z_m, \
				cellx_x_m, cellx_y_m, cellx_z_m, fccx_m, fccy_m, fccz_m, \
				rface0_m, rface1_m, facn_m, edgeNum, vertexNum);\
	swCalcLudsFcc_master(cOpt_m, massFlux_m, facex_x_m, facex_y_m, facex_z_m,\
			   	cellx_x_m, cellx_y_m, cellx_z_m, fccx_m, fccy_m, fccz_m, \
				rface0_m, rface1_m, facn_m, CS, edgeNum, vertexNum);\
	printf("test\n"); \
	getTime(time1);\
	for(int iOpt=0;iOpt<optNum;iOpt++) \
	{ \
		backEdgeData_m  = *(cOpt_m[iOpt].data_p->backEdgeData); \
		frontEdgeData_m = *(cOpt_m[iOpt].data_p->frontEdgeData); \
		selfConnData_m  = *(cOpt_m[iOpt].data_p->selfConnData); \
		vertexData_m    = *(cOpt_m[iOpt].data_p->vertexData); \
		cOpt_m[iOpt].fun_host(&backEdgeData_m, &frontEdgeData_m, \
					&selfConnData_m, &vertexData_m, &paraData, \
					rowAddr,colAddr, cOpt_m[iOpt].data);\
	} \
    getTime(time2);\
	printf("CPU Processor Time: %f us\n", (time2-time1)*1000000);\
	CG_init();\
    getTime(time1);\
	swCalcLudsFcc_host(&mlbIter, massFlux, facex_x, facex_y, facex_z, cellx_x, \
				cellx_y, cellx_z, fccx, fccy, fccz, rface0, rface1, facn, \
			   	CS, gamblend, edgeNum, vertexNum);\
	getTime(time2);\
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000);\
	CG_halt();\
\
	checkResult(fccx_m, fccx, edgeNum*dims);\
	checkResult(fccy_m, fccy, edgeNum*dims);\
	checkResult(fccz_m, fccz, edgeNum*dims);\
	checkResult(rface0_m, rface0, edgeNum*dims);\
	checkResult(rface1_m, rface1, edgeNum*dims);\
}

#ifdef __cplusplus
}
#endif

#endif
