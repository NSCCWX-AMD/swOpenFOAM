#ifndef INTERPOLATEVISCOSITY_H
#define INTERPOLATEVISCOSITY_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swInterpolateViscosity);
define_e2v_slaveFunPtr(swInterpolateViscosity);
define_e2v_FunPtr(swInterpolateViscosity_1);
define_e2v_slaveFunPtr(swInterpolateViscosity_1);

#define swInterpolateViscosity_data(visSx, visSy, visSz, visPNx, visPNy, \
			visPNz, rface0, rface1, face_nx, face_ny, face_nz, face_dx, \
			face_dy, face_dz, face_d, viseff, facn, edgeNum, vertexNum) \
{ \
	int i,j; \
	visSx  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visSy  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visSz  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visPNx = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visPNy = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visPNz = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface0 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface1 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* phi[8] = {visSx, visSy, visSz, visPNx, visPNy, visPNz, \
		rface0, rface1}; \
	for(i=0;i<8;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			phi[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	face_nx = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_ny = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_nz = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_dx = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_dy = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_dz = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_d  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	facn    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* S[8] = {face_nx, face_ny, face_nz, face_dx, face_dy, face_dz, \
		face_d, facn}; \
	for(i=0;i<8;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			S[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
	viseff = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	for(j=0;j<vertexNum;j++) \
	{ \
		viseff[j] = (swFloat)(j+1)/(j+2); \
	} \
} 

#define swInterpolateViscosity_master(cOpt, visSx, visSy, visSz, visPNx, \
			visPNy, visPNz, rface0, rface1, face_nx, face_ny, face_nz, \
			face_dx, face_dy, face_dz, face_d, viseff, facn, \
			edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	swFloat* phi[8] = {visSx, visSy, visSz, visPNx, visPNy, visPNz, \
		rface0, rface1}; \
	swFloat* S[8] = {face_nx, face_ny, face_nz, face_dx, face_dy, face_dz, \
		face_d, face_d}; \
	Arrays backEdgeData[8], frontEdgeData[8]; \
	Arrays selfConnData[8], vertexData[8]; \
	FieldData data[8]; \
	Arrays backEdgeData_p[8], frontEdgeData_p[8]; \
	Arrays selfConnData_p[8], vertexData_p[8]; \
	FieldData data_p[8]; \
	for(iArray=0;iArray<8;iArray++) \
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
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_nx);\
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visSx); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, facn);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, COPYIN, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity; \
	cOpt[opt].fun_host = swInterpolateViscosity; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_ny);\
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visSy); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity; \
	cOpt[opt].fun_host = swInterpolateViscosity; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_nz);\
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visSz); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity; \
	cOpt[opt].fun_host = swInterpolateViscosity; \
	\
	/* operator 3 */ \
	opt = 3; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_dx);\
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visPNx); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity; \
	cOpt[opt].fun_host = swInterpolateViscosity; \
	\
	/* operator 4 */ \
	opt = 4; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_dy);\
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visPNy); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity; \
	cOpt[opt].fun_host = swInterpolateViscosity; \
	\
	/* operator 5 */ \
	opt = 5; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_dz);\
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visPNz); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity; \
	cOpt[opt].fun_host = swInterpolateViscosity; \
	\
	/* operator 6 */ \
	opt = 6; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface0); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, face_d);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity_1; \
	cOpt[opt].fun_host = swInterpolateViscosity_1; \
	\
	/* operator 7 */ \
	opt = 7; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface1); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);\
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, face_d);\
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);\
	cOpt[opt].fun_slave = slave_swInterpolateViscosity_1; \
	cOpt[opt].fun_host = swInterpolateViscosity_1; \
	\
}

#define swInterpolateViscosity_test() \
{ \
    int optNum = 8; \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *visSx, *visSy, *visSz, *visPNx, *visPNy, *visPNz; \
	swFloat *rface0, *rface1, *face_nx, *face_ny, *face_nz; \
	swFloat *face_dx, *face_dy, *face_dz, *face_d, *viseff, *facn; \
	swFloat *visSx_m, *visSy_m, *visSz_m, *visPNx_m, *visPNy_m, *visPNz_m; \
	swFloat *rface0_m, *rface1_m, *face_nx_m, *face_ny_m, *face_nz_m; \
	swFloat *face_dx_m, *face_dy_m, *face_dz_m, *face_d_m; \
	swFloat *viseff_m, *facn_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	swInterpolateViscosity_data(visSx, visSy, visSz, visPNx, visPNy, \
				visPNz, rface0, rface1, face_nx, face_ny, face_nz, \
				face_dx, face_dy, face_dz, face_d, viseff, facn, \
				edgeNum, vertexNum); \
    swInterpolateViscosity_data(visSx_m, visSy_m, visSz_m, visPNx_m, \
				visPNy_m, visPNz_m, rface0_m, rface1_m, face_nx_m, \
				face_ny_m, face_nz_m, face_dx_m, face_dy_m, face_dz_m, \
				face_d_m, viseff_m, facn_m, edgeNum, vertexNum); \
	swInterpolateViscosity_master(cOpt_m, visSx_m, visSy_m, visSz_m, \
				visPNx_m, visPNy_m, visPNz_m, rface0_m, rface1_m, \
				face_nx_m, face_ny_m, face_nz_m, face_dx_m, face_dy_m, \
				face_dz_m, face_d_m, viseff_m, facn_m, \
				edgeNum, vertexNum); \
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
	swInterpolateViscosity_host(&mlbIter, visSx, visSy, visSz, visPNx, visPNy, \
				visPNz, rface0, rface1, face_nx, face_ny, face_nz, \
				face_dx, face_dy, face_dz, face_d, viseff, facn, \
				edgeNum, vertexNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(visSx_m, visSx, edgeNum*dims); \
	checkResult(visSy_m, visSy, edgeNum*dims); \
	checkResult(visSz_m, visSz, edgeNum*dims); \
	checkResult(visPNx_m, visPNx, edgeNum*dims); \
	checkResult(visPNy_m, visPNy, edgeNum*dims); \
	checkResult(visPNz_m, visPNz, edgeNum*dims); \
	checkResult(rface0_m, rface0, edgeNum*dims); \
	checkResult(rface1_m, rface1, edgeNum*dims); \
}

#ifdef __cplusplus
}
#endif

#endif
