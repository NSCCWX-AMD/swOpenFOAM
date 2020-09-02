#ifndef CALCULATEVISFLUX_H
#define CALCULATEVISFLUX_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swCalculateVisFlux);
define_e2v_slaveFunPtr(swCalculateVisFlux);
define_e2v_FunPtr(swCalculateVisFlux_multiply2);
define_e2v_slaveFunPtr(swCalculateVisFlux_multiply2);

#define swCalculateVisFlux_data(facn, visSy,dudy, dvdx, dvdy, dvdz,dwdy, \
			Su, Sv, Sw, dudx, dudz, dwdx, dwdz, visSx, visSz, visPNx, \
			visPNy, visPNz, edgeNum, vertexNum) \
{ \
	int i,j; \
	dudy = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dvdx = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dvdy = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dvdz = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dwdy = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dudx = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dudz = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dwdx = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dwdz = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* phi[9] = {dudy, dvdx, dvdy, dvdz, dwdy, dudx, dudz, dwdx, \
		dwdz}; \
	for(i=0;i<9;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			phi[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	Su = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	Sv = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	Sw = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* S[3] = {Su, Sv, Sw}; \
	for(i=0;i<3;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			S[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
	facn   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visSx  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visSy  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visSz  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visPNx = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visPNy = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visPNz = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat *edge[7] = {facn, visSx, visSy, visSz, visPNx, visPNy, visPNz};\
	for(i=0;i<7;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = (swFloat)(3*j+i+1)/(3*j+i+2); \
		} \
	} \
} 

#define swCalculateVisFlux_master(cOpt, facn, visSy, dudy, dvdx, dvdy, \
			dvdz, dwdy, Su, Sv, Sw, dudx, dudz, dwdx, dwdz, visSx, visSz, \
			visPNx, visPNy, visPNz, edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	Arrays backEdgeData[24], frontEdgeData[24]; \
	Arrays selfConnData[24], vertexData[24]; \
	FieldData data[24]; \
	Arrays backEdgeData_p[24], frontEdgeData_p[24]; \
	Arrays selfConnData_p[24], vertexData_p[24]; \
	FieldData data_p[24]; \
	for(iArray=0;iArray<24;iArray++) \
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
		\
	} \
	/* operator 0 */ \
	opt = 0; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visSy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux_multiply2; \
	cOpt[opt].fun_host  = swCalculateVisFlux_multiply2; \
	\
	/* operator 3 */ \
	opt = 3; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum,UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum,UPDATED, visSy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 4 */ \
	opt = 4; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 5 */ \
	opt = 5; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visSx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux_multiply2; \
	cOpt[opt].fun_host  = swCalculateVisFlux_multiply2; \
	\
	/* operator 6 */ \
	opt = 6; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 7 */ \
	opt = 7; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 8 */ \
	opt = 8; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 9 */ \
	opt = 9; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 10 */ \
	opt = 10; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visSz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 11 */ \
	opt = 11; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 12 */ \
	opt = 12; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 13 */ \
	opt = 13; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 14 */ \
	opt = 14; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux_multiply2; \
	cOpt[opt].fun_host  = swCalculateVisFlux_multiply2; \
	\
	/* operator 15 */ \
	opt = 15; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visPNx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 16 */ \
	opt = 16; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 17 */ \
	opt = 17; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 18 */ \
	opt = 18; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visPNy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 19 */ \
	opt = 19; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 20 */ \
	opt = 20; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 21 */ \
	opt = 21; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visPNz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 22 */ \
	opt = 22; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
	/* operator 23 */ \
	opt = 23; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdz); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateVisFlux; \
	cOpt[opt].fun_host  = swCalculateVisFlux; \
	\
}

#define swCalculateVisFlux_test() \
{ \
    int optNum = 24; \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
    swFloat *dudy, *dvdx, *dvdy, *dvdz, *dwdy, *Su, *Sv, *Sw; \
	swFloat *facn, *visSy, *dudx, *dudz, *dwdx, *dwdz, *visSx, *visSz; \
	swFloat *visPNx, *visPNy, *visPNz; \
	swFloat *dudy_m, *dvdx_m, *dvdy_m, *dvdz_m, *dwdy_m; \
	swFloat *Su_m, *Sv_m, *Sw_m, *facn_m, *visSy_m; \
	swFloat *dudx_m, *dudz_m, *dwdx_m, *dwdz_m, *visSx_m, *visSz_m; \
	swFloat *visPNx_m, *visPNy_m, *visPNz_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	swCalculateVisFlux_data(facn, visSy, dudy, dvdx, dvdy, dvdz, dwdy, \
				Su, Sv, Sw, dudx, dudz, dwdx, dwdz, visSx, visSz, visPNx, \
				visPNy, visPNz, edgeNum, vertexNum); \
	swCalculateVisFlux_data(facn_m, visSy_m, dudy_m, dvdx_m, dvdy_m, \
				dvdz_m, dwdy_m, Su_m, Sv_m, Sw_m, dudx_m, dudz_m, dwdx_m, \
				dwdz_m, visSx_m, visSz_m, visPNx_m, visPNy_m, visPNz_m, \
				edgeNum, vertexNum); \
	swCalculateVisFlux_master(cOpt_m, facn_m, visSy_m, dudy_m, dvdx_m, \
				dvdy_m, dvdz_m, dwdy_m, Su_m, Sv_m, Sw_m, dudx_m, dudz_m, \
				dwdx_m, dwdz_m, visSx_m, visSz_m, visPNx_m, visPNy_m, \
				visPNz_m, edgeNum, vertexNum); \
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
	swCalculateVisFlux_host(&mlbIter, facn, visSy, dudy, dvdx, dvdy, dvdz, \
				dwdy, Su, Sv, Sw, dudx, dudz, dwdx, dwdz, visSx, visSz, \
				visPNx, visPNy, visPNz, edgeNum, vertexNum); \
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
