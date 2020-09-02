#ifndef SEPARATE_VECTOR_H
#define SEPARATE_VECTOR_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_array_FunPtr(swSeparate_vector);
define_array_slaveFunPtr(swSeparate_vector);
//define_array_FunPtr(swSeparate_vector_2);
//define_array_slaveFunPtr(swSeparate_vector_2);
//define_array_FunPtr(swSeparate_UDS_thermal);
//define_array_slaveFunPtr(swSeparate_UDS_thermal);

#define swSeparate_vector_data(mf, phiUDS, face_x0, dPhidXac0, \
			face_x1, dPhidXac1, face_x2, dPhidXac2, XU0, XU1, XU2, visac, \
			face_D, Xpn0, Xpn1, Xpn2, face_n0, face_n1, face_n2, dPhidXU0, \
			dPhidXU1, dPhidXU2, fu, rface_1, rface_2, rcpac, denac, \
			Xac0, Xac1, Xac2, phiCDS, phiCDS2, edgeNum, vertexNum) \
{ \
	int i,j; \
	mf     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiUDS = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x0 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac0 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac1 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac2 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU0       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU1       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU2       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_D    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xpn0      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xpn1      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xpn2      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n0   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU0  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU1  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU2  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fu        = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface_1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface_2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rcpac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	denac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xac0      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xac1      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xac2      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiCDS    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiCDS2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* edge[32] = {mf, phiUDS, face_x0, dPhidXac0, face_x1, \
		dPhidXac1, face_x2, dPhidXac2, XU0, XU1, XU2, visac, face_D, \
		Xpn0, Xpn1, Xpn2, face_n0, face_n1, face_n2, dPhidXU0, dPhidXU1, \
		dPhidXU2, fu, rface_1, rface_2, rcpac, denac, Xac0, Xac1, \
		Xac2, phiCDS, phiCDS2} ;\
	for(i=0;i<32;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = (swFloat)(2*i+j+1)/(2*i+j+2); \
		} \
	} \
} 

#define swSeparate_vector_master(cOpt, mf, phiUDS, face_x0, dPhidXac0, \
			face_x1, dPhidXac1, face_x2, dPhidXac2, XU0, XU1, XU2, visac, \
			face_D, Xpn0, Xpn1, Xpn2, face_n0, face_n1, face_n2, dPhidXU0, \
			dPhidXU1, dPhidXU2, fu, rface_1, rface_2, rcpac, denac, \
			Xac0, Xac1, Xac2, phiCDS, phiCDS2, phiCDDedgeNum, vertexNum) \
{ \
	int iArray; \
	Arrays backEdgeData, frontEdgeData; \
	Arrays selfConnData, vertexData; \
	FieldData data; \
	constructSingleArray(frontEdgeData, 1, edgeNum, COPYOUT, fu); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_x0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, XU0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXU0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_x1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, XU1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXU1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_x2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, XU2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXU2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, mf); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, phiUDS); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, visac); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_D); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xpn0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXac0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xpn1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXac1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xpn2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXac2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_n0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_n1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_n2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYOUT, rface_1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYOUT, rface_2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, rcpac); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, denac); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xac0); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xac1); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xac2); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, phiCDS); \
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, phiCDS2); \
	constructEmptyArray(backEdgeData); \
	constructEmptyArray(selfConnData); \
	constructEmptyArray(vertexData); \
	data.backEdgeData  = &backEdgeData; \
	data.frontEdgeData = &frontEdgeData; \
	data.selfConnData  = &selfConnData; \
	data.vertexData    = &vertexData; \
	cOpt[0].fun_slave = slave_swSeparate_vector; \
	cOpt[0].fun_host  = swSeparate_vector; \
	cOpt[0].data = &data; \
}

#define swSeparate_vector_test() \
{ \
    int optNum = 1; \
	coupledArrayOperator *cOpt_m \
		=(coupledArrayOperator*)malloc(optNum*sizeof(coupledArrayOperator)); \
	swFloat *mf, *phiUDS, *face_x0, *Xac0, *dPhidXac0, *face_x1, *Xac1; \
	swFloat *dPhidXac1, *face_x2, *Xac2, *dPhidXac2, *phiCDS, *phiCDS2; \
	swFloat *XU0, *XU1, *XU2, *rcpac, *visac, *denac, *face_D, *Xpn0; \
	swFloat *Xpn1, *Xpn2, *face_n0, *face_n1, *face_n2, *fu; \
	swFloat *dPhidXU0, *dPhidXU1, *dPhidXU2, *rface_1, *rface_2; \
	swFloat *mf_m, *phiUDS_m, *face_x0_m, *Xac0_m, *dPhidXac0_m, \
		*face_x1_m, *Xac1_m; \
	swFloat *dPhidXac1_m, *face_x2_m, *Xac2_m, *dPhidXac2_m, *phiCDS_m, \
		*phiCDS2_m; \
	swFloat *XU0_m, *XU1_m, *XU2_m, *rcpac_m, *visac_m, *denac_m, \
		*face_D_m, *Xpn0_m; \
	swFloat *Xpn1_m, *Xpn2_m, *face_n0_m, *face_n1_m, *face_n2_m, *fu_m; \
	swFloat *dPhidXU0_m, *dPhidXU1_m, *dPhidXU2_m, *rface_1_m, *rface_2_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	swFloat cvScheme       = 2; \
	swFloat vis_type       = 2; \
	swFloat csBlend        = 1.1111111; \
	swFloat vis_lam        = 3.3333333; \
	swFloat vis_sigma      = 4.4444444; \
	swFloat vis_PrScNr     = 6.6666666; \
	swFloat vis_lambda     = 7.4354354; \
	swFloat vis_diffuse    = 6.6662342; \
	swFloat csFaceCorrect  = 4.5662342; \
	swFloat paras[9] = {csBlend, vis_lam, vis_sigma, vis_PrScNr, \
				vis_lambda, vis_diffuse, vis_type, cvScheme, \
				csFaceCorrect};\
	constructSingleArray(paraData,1,9,COPYIN,&paras[0]); \
	swSeparate_vector_data(mf, phiUDS, face_x0, dPhidXac0, face_x1, \
				dPhidXac1, face_x2, dPhidXac2, XU0, XU1, XU2, visac, \
				face_D, Xpn0, Xpn1, Xpn2, face_n0, face_n1, face_n2, \
				dPhidXU0, dPhidXU1,	dPhidXU2, fu, rface_1, rface_2, \
				rcpac, denac, Xac0, Xac1, Xac2, phiCDS, phiCDS2, \
				edgeNum, vertexNum); \
	swSeparate_vector_data(mf_m, phiUDS_m, face_x0_m, dPhidXac0_m, \
				face_x1_m, dPhidXac1_m, face_x2_m, dPhidXac2_m, XU0_m, \
				XU1_m, XU2_m, visac_m, face_D_m, Xpn0_m, Xpn1_m, Xpn2_m,\
			   	face_n0_m, face_n1_m, face_n2_m, dPhidXU0, dPhidXU1, \
				dPhidXU2, fu_m, rface_1_m, rface_2_m, rcpac_m, denac_m, \
				Xac0_m, Xac1_m, Xac2_m, phiCDS_m, phiCDS2_m, \
				edgeNum, vertexNum); \
	swSeparate_vector_master(cOpt_m, mf_m, phiUDS_m, face_x0_m, dPhidXac0_m, \
				face_x1_m, dPhidXac1_m, face_x2_m, dPhidXac2_m, XU0_m, \
				XU1_m, XU2_m, visac_m, face_D_m, Xpn0_m, Xpn1_m, Xpn2_m, \
				face_n0_m, face_n1_m, face_n2_m, dPhidXU0, dPhidXU1, \
				dPhidXU2, fu_m, rface_1_m, rface_2_m, rcpac_m, denac_m, \
				Xac0_m, Xac1_m, Xac2_m, phiCDS_m, phiCDS2_m, \
				edgeNum, vertexNum); \
	printf("test\n"); \
	getTime(time1); \
	for(int iOpt=0;iOpt<optNum;iOpt++) \
	{ \
		cOpt_m[iOpt].fun_host(&backEdgeData_m, &frontEdgeData_m, \
					&selfConnData_m, &vertexData_m, &paraData, \
					cOpt_m[iOpt].data); \
	} \
    getTime(time2); \
	printf("CPU Processor Time: %f us\n", (time2-time1)*1000000); \
	CG_init(); \
    getTime(time1); \
	swSeparate_vector_host(&mlbIter, mf, phiUDS, face_x0, dPhidXac0, face_x1, \
				dPhidXac1, face_x2, dPhidXac2, XU0, XU1, XU2, visac, \
			   	face_D, Xpn0, Xpn1, Xpn2, face_n0, face_n1, face_n2, \
				dPhidXU0, dPhidXU1, dPhidXU2, fu, rface_1, rface_2, \
				rcpac, denac, Xac0, Xac1, Xac2, phiCDS, phiCDS2, csBlend,\
				vis_lam, vis_sigma, vis_PrScNr, vis_lambda, vis_diffuse, \
				vis_type, cvScheme, csFaceCorrect, edgeNum, vertexNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(fu_m, fu, edgeNum); \
	checkResult(rface_1_m, rface_1, edgeNum); \
	checkResult(rface_2_m, rface_2, edgeNum); \
}

#ifdef __cplusplus
}
#endif

#endif
