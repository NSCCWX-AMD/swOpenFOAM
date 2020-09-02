#ifndef CALCULATEFCCFLUX_H
#define CALCULATEFCCFLUX_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swCalculateFccFlux);
define_e2v_slaveFunPtr(swCalculateFccFlux);

#define swCalculateFccFlux_data(massFlux, fccx, dudx, dvdx, dwdx, \
			Su, Sv, Sw, fccy, fccz, dudy, dvdy, dwdy, dudz, dvdz, dwdz, \
			edgeNum, vertexNum) \
{ \
	int i,j; \
	dudx = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dvdx = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dwdx = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dudy = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dvdy = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dwdy = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dudz = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dvdz = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	dwdz = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* phi[9] = {dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, \
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
	massFlux  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fccx      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fccy      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fccz      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	for(j=0;j<edgeNum;j++) \
	{ \
		massFlux[j] = (swFloat)(j+1)/(j+2); \
	} \
	for(j=0;j<edgeNum;j++) \
	{ \
		fccx[j] = (swFloat)(j+1)/(j+2); \
		fccy[j] = (swFloat)(j+2)/(j+3); \
		fccz[j] = (swFloat)(j+3)/(j+2); \
	} \
} 

#define swCalculateFccFlux_master(cOpt, massFlux, fccx, dudx, dvdx, dwdx, \
			Su, Sv, Sw, fccy, fccz, dudy, dvdy, dwdy, dudz, dvdz, dwdz, \
			edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	Arrays backEdgeData[9], frontEdgeData[9]; \
	Arrays selfConnData[9], vertexData[9]; \
	FieldData data[9]; \
	Arrays backEdgeData_p[9], frontEdgeData_p[9]; \
	Arrays selfConnData_p[9], vertexData_p[9]; \
	FieldData data_p[9]; \
	for(iArray=0;iArray<9;iArray++) \
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
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudx); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, COPYIN, fccx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccx); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 3 */ \
	opt = 3; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, COPYIN, fccy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 4 */ \
	opt = 4; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdy); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 5 */ \
	opt = 5; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccy); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 6 */ \
	opt = 6; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, COPYIN, fccz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 7 */ \
	opt = 7; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sv); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
	/* operator 8 */ \
	opt = 8; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdz); \
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sw); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux); \
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccz); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swCalculateFccFlux; \
	cOpt[opt].fun_host  = swCalculateFccFlux; \
	\
}

#define swCalculateFccFlux_test() \
{ \
    int optNum = 9; \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
    swFloat *dudx, *dvdx, *dwdx, *Su, *Sv, *Sw, *dudy, *dvdy, *dwdy; \
	swFloat *massFlux, *fccx, *fccy, *fccz, *dudz, *dvdz, *dwdz; \
	swFloat *dudx_m, *dvdx_m, *dwdx_m, *dudy_m, *dvdy_m, *dwdy_m; \
	swFloat *Su_m, *Sv_m, *Sw_m, *massFlux_m, *fccx_m, *fccy_m, *fccz_m; \
	swFloat *dudz_m, *dvdz_m, *dwdz_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	swCalculateFccFlux_data(massFlux, fccx, dudx, dvdx, dwdx, \
				Su, Sv, Sw, fccy, fccz, dudy, dvdy, dwdy, dudz, dvdz, dwdz,\
				edgeNum, vertexNum); \
	swCalculateFccFlux_data(massFlux_m, fccx_m, dudx_m, dvdx_m, dwdx_m, \
				Su_m, Sv_m, Sw_m, fccy_m, fccz_m, dudy_m, dvdy_m, dwdy_m, \
				dudz_m, dvdz_m, dwdz_m, edgeNum, vertexNum); \
	swCalculateFccFlux_master(cOpt_m, massFlux_m, fccx_m, dudx_m, dvdx_m, \
				dwdx_m, Su_m, Sv_m, Sw_m, fccy_m, fccz_m, dudy_m, dvdy_m, \
				dwdy_m, dudz_m, dvdz_m, dwdz_m, edgeNum, vertexNum); \
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
	swCalculateFccFlux_host(&mlbIter, massFlux, fccx, dudx, dvdx, dwdx, \
				Su, Sv, Sw, fccy, fccz, dudy, dvdy, dwdy, dudz, dvdz, dwdz,\
				edgeNum, vertexNum); \
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
