#ifndef MERGE4_INTERPOLATION_H
#define MERGE4_INTERPOLATION_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swMerge4_lam);
define_e2v_slaveFunPtr(swMerge4_lam);
define_e2v_FunPtr(swMerge4_mf);
define_e2v_slaveFunPtr(swMerge4_mf);
define_e2v_FunPtr(swMerge4_dPhi);
define_e2v_slaveFunPtr(swMerge4_dPhi);

#define swMerge4_intp_data(face_lam, vis, visac, gradPhi0, dPhidXac0, \
		    gradPhi1, dPhidXac1, gradPhi2, dPhidXac2, mf, cell_x0, XU0, \
			cell_x1, XU1, cell_x2, XU2, dPhidXU0, dPhidXU1, dPhidXU2, \
			phi, phiUDS, dPhi, edgeNum, vertexNum) \
{ \
	int i,j; \
	face_lam  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac0 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac1 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXac2 = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	mf        = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU0       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU1       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XU2       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU0  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU1  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidXU2  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiUDS    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhi      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* edge[14] = {face_lam, visac, dPhidXac0, dPhidXac1, dPhidXac2, \
		mf, XU0, XU1, XU2, dPhidXU0, dPhidXU1, dPhidXU2, phiUDS, dPhi}; \
	for(i=0;i<14;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	vis      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	gradPhi0 = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	gradPhi1 = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	gradPhi2 = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cell_x0  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cell_x1  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cell_x2  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	phi      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* cell[8] = {vis, gradPhi0, gradPhi1, gradPhi2, cell_x0, \
		cell_x1, cell_x2, phi}; \
	for(i=0;i<8;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			cell[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
} 

#define swMerge4_intp_host(cOpt, face_lam, vis, visac, gradPhi0, dPhidXac0,\
		    gradPhi1, dPhidXac1, gradPhi2, dPhidXac2, mf, cell_x0, XU0, \
			cell_x1, XU1, cell_x2, XU2, dPhidXU0, dPhidXU1, dPhidXU2, \
			phi, phiUDS, dPhi, edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	Arrays backEdgeData[14], frontEdgeData[14]; \
	Arrays selfConnData[14], vertexData[14]; \
	FieldData data[14]; \
	Arrays backEdgeData_p[14], frontEdgeData_p[14]; \
	Arrays selfConnData_p[14], vertexData_p[14]; \
	FieldData data_p[14]; \
	for(iArray=0;iArray<14;iArray++) \
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
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visac); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, vis); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_lam; \
	cOpt[opt].fun_host  = swMerge4_lam; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXac0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_lam; \
	cOpt[opt].fun_host  = swMerge4_lam; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXac1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_lam; \
	cOpt[opt].fun_host  = swMerge4_lam; \
	\
	/* operator 3 */ \
	opt = 3; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXac2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_lam; \
	cOpt[opt].fun_host  = swMerge4_lam; \
	\
	/* operator 4 */ \
	opt = 4; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 5 */ \
	opt = 5; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 6 */ \
	opt = 6; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 7 */ \
	opt = 7; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 8 */ \
	opt = 8; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 9 */ \
	opt = 9; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 10 */ \
	opt = 10; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				phiUDS); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_mf; \
	cOpt[opt].fun_host  = swMerge4_mf; \
	\
	/* operator 11 */ \
	opt = 11; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, \
				dPhidXac0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, dPhi); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_dPhi; \
	cOpt[opt].fun_host  = swMerge4_dPhi; \
	\
	/* operator 12 */ \
	opt = 12; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, \
				dPhidXac1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, dPhi); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_dPhi; \
	cOpt[opt].fun_host  = swMerge4_dPhi; \
	\
	/* operator 13 */ \
	opt = 13; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, \
				dPhidXac2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYOUT, dPhi); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge4_dPhi; \
	cOpt[opt].fun_host  = swMerge4_dPhi; \
	\
}

#define swMerge4_intp_test() \
{ \
    int optNum = 14; \
	coupledOperator *cOpt \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *face_lam, *vis, *visac, *gradPhi0, *dPhidXac0, *gradPhi1; \
	swFloat *dPhidXac1, *gradPhi2, *dPhidXac2; \
	swFloat *mf, *cell_x0, *XU0, *cell_x1, *XU1, *cell_x2, *XU2; \
	swFloat *dPhidXU0, *dPhidXU1, *dPhidXU2, *phi, *phiUDS, *dPhi; \
	swFloat *face_lam_m, *vis_m, *visac_m, *gradPhi0_m, *dPhidXac0_m; \
	swFloat *gradPhi1_m, *dPhidXac1_m, *gradPhi2_m, *dPhidXac2_m; \
	swFloat *mf_m, *cell_x0_m, *XU0_m, *cell_x1_m, *XU1_m, *cell_x2_m; \
   	swFloat *XU2_m, *dPhi_m; \
	swFloat *dPhidXU0_m, *dPhidXU1_m, *dPhidXU2_m, *phi_m, *phiUDS_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	swMerge4_intp_data(face_lam, vis, visac, gradPhi0, dPhidXac0, gradPhi1,\
				dPhidXac1, gradPhi2, dPhidXac2, mf, cell_x0, XU0, cell_x1, \
				XU1, cell_x2, XU2, dPhidXU0, dPhidXU1, dPhidXU2, phi, \
				phiUDS, dPhi, edgeNum, vertexNum); \
	swMerge4_intp_host(cOpt, face_lam, vis, visac, gradPhi0, dPhidXac0, \
			   	gradPhi1, dPhidXac1, gradPhi2, dPhidXac2, mf, cell_x0, \
				XU0, cell_x1, XU1, cell_x2, XU2,dPhidXU0, dPhidXU1, \
			   	dPhidXU2, phi, phiUDS, dPhi, edgeNum, vertexNum); \
	swMerge4_intp_data(face_lam_m, vis_m, visac_m, gradPhi0_m, dPhidXac0_m,\
				gradPhi1_m, dPhidXac1_m, gradPhi2_m, dPhidXac2_m, mf_m, \
				cell_x0_m, XU0_m, cell_x1_m, XU1_m, cell_x2_m, XU2_m, \
				dPhidXU0_m, dPhidXU1_m, dPhidXU2_m, phi_m, phiUDS_m, \
				dPhi_m, edgeNum, vertexNum); \
	swMerge4_intp_host(cOpt_m, face_lam_m, vis_m, visac_m, gradPhi0_m, \
				dPhidXac0_m, gradPhi1_m, dPhidXac1_m, gradPhi2_m, \
				dPhidXac2_m, mf_m, cell_x0_m, XU0_m, cell_x1_m, XU1_m, \
				cell_x2_m, XU2_m, dPhidXU0_m, dPhidXU1_m, dPhidXU2_m, \
				phi_m, phiUDS_m, dPhi_m, edgeNum, vertexNum); \
	optNum =14; \
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
 	mlbIter.edge2VertexIteration(&paraData, cOpt, optNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(visac_m, visac, edgeNum); \
	checkResult(dPhidXac0_m, dPhidXac0, edgeNum); \
	checkResult(dPhidXac1_m, dPhidXac1, edgeNum); \
	checkResult(dPhidXac2_m, dPhidXac2, edgeNum); \
	checkResult(XU0_m, XU0, edgeNum); \
	checkResult(XU1_m, XU1, edgeNum); \
	checkResult(XU2_m, XU2, edgeNum); \
	checkResult(dPhidXU0_m, dPhidXU0, edgeNum); \
	checkResult(dPhidXU1_m, dPhidXU1, edgeNum); \
	checkResult(dPhidXU2_m, dPhidXU2, edgeNum); \
	checkResult(dPhi_m, dPhi, edgeNum); \
}

#ifdef __cplusplus
}
#endif

#endif
