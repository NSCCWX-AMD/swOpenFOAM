#ifndef UMBT_CONV_H
#define UMBT_CONV_H
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
//define 2 function pointers
define_e2v_FunPtr(umbt_conv);
define_e2v_slaveFunPtr(umbt_conv);

#define umbt_conv_data(n_vg, u, cen, gra, rhs, edgeNum, vertexNum) \
{ \
	int i,j,idim; \
	n_vg = (swFloat*)malloc(edgeNum*8*sizeof(swFloat)); \
	cen  = (swFloat*)malloc(edgeNum*3*sizeof(swFloat)); \
	for(i=0;i<edgeNum;i++) \
	{ \
		for(idim=0;idim<8;idim++) \
		{ \
			n_vg[i*8+idim] = (i+idim+2)/(i+idim+1); \
		} \
	} \
	for(i=0;i<edgeNum;i++) \
	{ \
		for(idim=0;idim<3;idim++) \
		{ \
			cen[i*3+idim] = (i+idim+1)/(i+idim+2); \
		} \
	} \
	\
	u   = (swFloat*)malloc(vertexNum*5*sizeof(swFloat)); \
	rhs = (swFloat*)malloc(vertexNum*5*sizeof(swFloat)); \
	gra = (swFloat*)malloc(vertexNum*19*sizeof(swFloat)); \
	for(i=0;i<vertexNum;i++) \
	{ \
		for(idim=0;idim<5;idim++) \
		{ \
			u[i*5+idim] = (2*i+idim+2)/(2*i+1); \
			rhs[i*5+idim] = (2*i+idim+1)/(2*i+2); \
		} \
	} \
	for(i=0;i<vertexNum;i++) \
	{ \
		for(idim=0;idim<19;idim++) \
		{ \
			gra[i*19+idim] = (2*i+idim+1)/(2*i+3); \
		} \
	} \
}

#define umbt_conv_master(cOpt, n_vg, u, cen, gra, rhs, edgeNum ,vertexNum) \
{ \
	Arrays backEdgeData,frontEdgeData,selfConnData,vertexData; \
	constructEmptyArray(backEdgeData); \
	constructSingleArray(frontEdgeData, 8, (edgeNum) ,COPYIN, n_vg); \
	addSingleArray(frontEdgeData, 3, (edgeNum), COPYIN, cen); \
	constructEmptyArray(selfConnData); \
	constructSingleArray(vertexData, 5, (vertexNum),COPYOUT,rhs); \
	addSingleArray(vertexData, 5, (vertexNum), COPYIN, u); \
	addSingleArray(vertexData, 19, (vertexNum), COPYIN, gra); \
	FieldData data \
		= {&backEdgeData, &frontEdgeData, &selfConnData, &vertexData}; \
	cOpt[0].fun_slave = slave_umbt_conv; \
	cOpt[0].fun_host  = umbt_conv; \
	cOpt[0].data      = &data; \
}

#define umbt_conv_test() \
{ \
	int optNum =1; \
	coupledOperator *cOpt = (coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	coupledOperator *cOpt_m = (coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *n_vg, *cen, *u, *rhs, *gra; \
	swFloat *n_vg_m, *cen_m, *u_m, *rhs_m, *gra_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	umbt_conv_data(n_vg, u, cen, gra, rhs, edgeNum, vertexNum); \
	umbt_conv_master(cOpt, n_vg, u, cen, gra, rhs, edgeNum, vertexNum); \
	umbt_conv_data(n_vg_m, u_m, cen_m, gra_m, rhs_m, edgeNum, vertexNum); \
	umbt_conv_master(cOpt_m, n_vg_m, u_m, cen_m, gra_m, rhs_m, \
			   	edgeNum, vertexNum); \
	printf("test\n"); \
    getTime(time1); \
	for(int iOpt=0;iOpt<optNum;iOpt++) \
	{ \
		cOpt_m[iOpt].fun_host(NULL,NULL,NULL,NULL,NULL, \
					startVertices,endVertices, cOpt_m[iOpt].data); \
	} \
	getTime(time2); \
	printf("CPU Processor Time: %f us\n", (time2-time1)*1000000); \
	CG_init(); \
	getTime(time1); \
	umbt_conv_host(&iter, n_vg, u, cen, gra, rhs, edgeNum, vertexNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(rhs_m, rhs, vertexNum*5); \
}

#ifdef __cplusplus
}
#endif

#endif
