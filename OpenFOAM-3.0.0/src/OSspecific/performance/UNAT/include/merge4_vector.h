#ifndef MERGE4_VECTOR_H
#define MERGE4_VECTOR_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_array_FunPtr(swMerge4_vector_1);
define_array_slaveFunPtr(swMerge4_vector_1);
define_array_FunPtr(swMerge4_vector_2);
define_array_slaveFunPtr(swMerge4_vector_2);
define_array_FunPtr(swMerge4_vector_3);
define_array_slaveFunPtr(swMerge4_vector_3);
define_array_FunPtr(swMerge4_vector_4);
define_array_slaveFunPtr(swMerge4_vector_4);
define_array_FunPtr(swMerge4_vector_5);
define_array_slaveFunPtr(swMerge4_vector_5);

#define swMerge4_vector_data(phiDelta, face_x0, XU0, dPhidXU0, face_x1, \
			XU1, dPhidXU1, face_x2, XU2, dPhidXU2, fde, visac, face_n0, \
			dPhidXac0, face_n1, dPhidXac1, face_n2, dPhidXac2, fdi, face_D,\
			dPhi, fce, mf, phiUDS, fci, edgeNum, vertexNum) \
{ \
	int i,j; \
	face_x0 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiDelta  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU0       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU1       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU2       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU0  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU1  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU2  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fde       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n0   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac0 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac1 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac2 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fdi       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_D    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhi      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fce       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	mf        = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiUDS    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fci       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* edge[25] = {face_x0, face_x1, face_x2, phiDelta, XU0, XU1, \
		XU2, dPhidXU0, dPhidXU1, dPhidXU2, fde, visac, face_n0, dPhidXac0,\
		dPhidXac1, dPhidXac2, face_n1, face_n2, fdi, face_D, dPhi, fce, \
		mf, phiUDS, fci} ;\
	for(i=0;i<25;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = (swFloat)(2*i+j+1)/(2*i+j+2); \
		} \
	} \
} 

#define swMerge4_vector_host(cOpt, phiDelta, face_x0, XU0, dPhidXU0, \
			face_x1, XU1, dPhidXU1, face_x2, XU2, dPhidXU2, fde, visac, \
			face_n0, dPhidXac0, face_n1, dPhidXac1, face_n2, dPhidXac2, \
			fdi, face_D, dPhi, fce, mf, phiUDS, fci, \
		   	edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	Arrays backEdgeData[9], frontEdgeData[9]; \
	Arrays selfConnData[9], vertexData[9]; \
	Arrays backEdgeData_p[9], frontEdgeData_p[9]; \
	Arrays selfConnData_p[9], vertexData_p[9]; \
	FieldData data[9]; \
	for(iArray=0;iArray<9;iArray++) \
	{ \
		constructEmptyArray(backEdgeData[iArray]); \
		constructEmptyArray(selfConnData[iArray]); \
		constructEmptyArray(vertexData[iArray]); \
		data[iArray].backEdgeData  = &backEdgeData[iArray]; \
		data[iArray].frontEdgeData = &frontEdgeData[iArray]; \
		data[iArray].selfConnData  = &selfConnData[iArray]; \
		data[iArray].vertexData    = &vertexData[iArray]; \
		cOpt[iArray].data = &data[iArray]; \
	} \
	opt =0; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, phiDelta); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_x0); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, XU0); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidXU0); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_1; \
	cOpt[opt].fun_host  = swMerge4_vector_1; \
	\
	opt =1; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, phiDelta); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_x1); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, XU1); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidXU1); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_1; \
	cOpt[opt].fun_host  = swMerge4_vector_1; \
	\
	opt =2; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, phiDelta); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_x2); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, XU2); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidXU2); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_1; \
	cOpt[opt].fun_host  = swMerge4_vector_1; \
	\
	opt =3; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, fde); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYINOUT, visac); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_n0); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidXac0); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_2; \
	cOpt[opt].fun_host  = swMerge4_vector_2; \
	\
	opt =4; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, fde); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, visac); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_n1); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidXac1); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_3; \
	cOpt[opt].fun_host  = swMerge4_vector_3; \
	\
	opt =5; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, fde); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, visac); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_n2); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidXac2); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_3; \
	cOpt[opt].fun_host  = swMerge4_vector_3; \
	\
	opt =6; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, fdi); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, visac); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_D); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhi); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_3; \
	cOpt[opt].fun_host  = swMerge4_vector_3; \
	\
	opt =7; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, fce); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, mf); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, phiUDS); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, phiDelta); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_4; \
	cOpt[opt].fun_host  = swMerge4_vector_4; \
	\
	opt =8; \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, fci); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, mf); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, phiUDS); \
	cOpt[opt].fun_slave = slave_swMerge4_vector_5; \
	cOpt[opt].fun_host  = swMerge4_vector_5; \
	\
}

#define swMerge4_vector_test() \
{ \
    int optNum = 9; \
	coupledArrayOperator *cOpt \
		=(coupledArrayOperator*)malloc(optNum*sizeof(coupledArrayOperator)); \
	coupledArrayOperator *cOpt_m \
		=(coupledArrayOperator*)malloc(optNum*sizeof(coupledArrayOperator)); \
	swFloat *phiDelta, *face_x0, *XU0, *dPhidXU0, *face_x1, *XU1; \
	swFloat *dPhidXU1, *face_x2, *XU2, *dPhidXU2, *fde, *visac, *face_n0; \
	swFloat *dPhidXac0, *face_n1, *dPhidXac1, *face_n2, *dPhidXac2; \
	swFloat *fdi, *face_D, *dPhi, *fce, *mf, *phiUDS, *fci; \
	swFloat *phiDelta_m, *face_x0_m, *XU0_m, *dPhidXU0_m, *face_x1_m; \
	swFloat *dPhidXU1_m, *face_x2_m, *XU2_m, *dPhidXU2_m, *XU1_m; \
	swFloat *fde_m, *visac_m, *face_n0_m, *dPhidXac0_m, *face_n1_m; \
	swFloat *face_n2_m, *dPhidXac2_m, *dPhidXac1_m, *fci_m; \
	swFloat *fdi_m, *face_D_m, *dPhi_m, *fce_m, *mf_m, *phiUDS_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	swFloat csBlend        = 1.1111111; \
	swFloat vis_lambda     = 2.2222222; \
	swFloat vis_lam        = 3.3333333; \
	swFloat vis_sigma      = 4.4444444; \
	swFloat csFaceCorrect  = 5.5555555; \
	swFloat vis_PrScNr     = 6.6666666; \
	swFloat vis_diffuse    = 7.7777777; \
	swFloat paras[7] = {csFaceCorrect, csBlend, vis_lambda, vis_lam, \
		vis_sigma, vis_diffuse, vis_PrScNr}; \
	constructSingleArray(paraData,1,7,COPYIN,&paras[0]); \
	swMerge4_vector_data(phiDelta, face_x0, XU0, dPhidXU0, face_x1, XU1, \
			   	dPhidXU1, face_x2, XU2, dPhidXU2, fde, visac, face_n0, \
				dPhidXac0, face_n1, dPhidXac1, face_n2, dPhidXac2, \
				fdi, face_D, dPhi, fce, mf, phiUDS, fci, \
				edgeNum, vertexNum); \
	swMerge4_vector_host(cOpt, phiDelta, face_x0, XU0, dPhidXU0, face_x1, \
				XU1, dPhidXU1, face_x2, XU2, dPhidXU2, fde, visac, face_n0,\
				dPhidXac0, face_n1, dPhidXac1, face_n2, dPhidXac2, \
				fdi, face_D, dPhi, fce, mf, phiUDS, fci, \
				edgeNum, vertexNum);\
	swMerge4_vector_data(phiDelta_m, face_x0_m, XU0_m, dPhidXU0_m, \
				face_x1_m, XU1_m, dPhidXU1_m, face_x2_m, XU2_m, dPhidXU2_m,\
			   	fde_m, visac_m, face_n0_m, dPhidXac0_m, face_n1_m, \
				dPhidXac1_m, face_n2_m, dPhidXac2_m, fdi_m, face_D_m, \
				dPhi_m, fce_m, mf_m, phiUDS_m, fci_m, \
				edgeNum, vertexNum); \
	swMerge4_vector_host(cOpt_m, phiDelta_m, face_x0_m, XU0_m, dPhidXU0_m, \
				face_x1_m, XU1_m, dPhidXU1_m, face_x2_m, XU2_m, dPhidXU2_m,\
				fde_m, visac_m, face_n0_m, dPhidXac0_m, face_n1_m, \
				dPhidXac1_m, face_n2_m, dPhidXac2_m, fdi_m, face_D_m, \
				dPhi_m, fce_m, mf_m, phiUDS_m, fci_m, \
				edgeNum, vertexNum); \
	optNum = 9; \
    for(int iOpt=0;iOpt<optNum;iOpt++) \
	{ \
		mlbIter.reorderEdgeData(cOpt[iOpt].data->backEdgeData, \
					cOpt[iOpt].data->frontEdgeData); \
		mlbIter.reorderVertexData(cOpt[iOpt].data->vertexData); \
		mlbIter.reorderVertexData(cOpt[iOpt].data->selfConnData); \
		mlbIter.reorderEdgeData(cOpt_m[iOpt].data->backEdgeData, \
					cOpt_m[iOpt].data->frontEdgeData); \
		mlbIter.reorderVertexData(cOpt_m[iOpt].data->vertexData); \
		mlbIter.reorderVertexData(cOpt_m[iOpt].data->selfConnData); \
	} \
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
 	mlbIter.arrayIteration(&paraData, cOpt, optNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(phiDelta_m, phiDelta, edgeNum); \
	checkResult(fde_m, fde, edgeNum); \
	checkResult(fdi_m, fdi, edgeNum); \
	checkResult(fce_m, fce, edgeNum); \
	checkResult(fci_m, fci, edgeNum); \
	checkResult(visac_m, visac, edgeNum); \
}

#ifdef __cplusplus
}
#endif

#endif
