#ifndef SEPARATE_INTERPOLATION_H
#define SEPARATE_INTERPOLATION_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swInterpolation_directed);
define_e2v_slaveFunPtr(swInterpolation_directed);
define_e2v_FunPtr(swInterpolation_weighted);
define_e2v_slaveFunPtr(swInterpolation_weighted);
define_e2v_FunPtr(swInterpolation_weightedSwap);
define_e2v_slaveFunPtr(swInterpolation_weightedSwap);
define_e2v_FunPtr(swInterpolation_constant_1);
define_e2v_slaveFunPtr(swInterpolation_constant_1);
define_e2v_FunPtr(swInterpolation_constant_0_5);
define_e2v_slaveFunPtr(swInterpolation_constant_0_5);

#define swSeparate_interpolation_data(face_lam, vis, visac, gradPhi0, \
			dPhidXac0, gradPhi1, dPhidXac1, gradPhi2, dPhidXac2, mf, \
			cell_x0, XU0, cell_x1, XU1, cell_x2, XU2, dPhidXU0, dPhidXU1, \
			dPhidXU2, phi, phiUDS, Xpn0, Xpn1, Xpn2, phiCDS, Xac0, \
			Xac1, Xac2, rcp, rcpac, den, denac, phiCDS2, \
			edgeNum, vertexNum) \
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
	Xpn0      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xpn1      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xpn2      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiCDS    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiCDS2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xac0      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xac1      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	Xac2      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rcpac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	denac     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* edge[23] = {face_lam, visac, dPhidXac0, dPhidXac1, dPhidXac2, \
		mf, XU0, XU1, XU2, dPhidXU0, dPhidXU1, dPhidXU2, phiUDS, \
		Xpn0, Xpn1, Xpn2, phiCDS, Xac0, Xac1, Xac2, rcpac, denac, phiCDS2};\
	for(i=0;i<23;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = j%(10+i); \
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
	rcp      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	den      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat* cell[10] = {vis, gradPhi0, gradPhi1, gradPhi2, cell_x0, \
		cell_x1, cell_x2, phi, rcp, den}; \
	for(i=0;i<10;i++) \
	{ \
		for(j=0;j<vertexNum;j++) \
		{ \
			cell[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
} 

#define swSeparate_interpolation_master(cOpt, face_lam, vis, visac, gradPhi0,\
		   	dPhidXac0, gradPhi1, dPhidXac1, gradPhi2, dPhidXac2, mf, \
			cell_x0, XU0, cell_x1, XU1, cell_x2, XU2, dPhidXU0, dPhidXU1, \
			dPhidXU2, phi, phiUDS, Xpn0, Xpn1, Xpn2, phiCDS, Xac0, Xac1, \
			Xac2, rcp, rcpac, den, denac, phiCDS2, edgeNum, vertexNum) \
{ \
	int iArray,opt; \
	Arrays backEdgeData[21], frontEdgeData[21]; \
	Arrays selfConnData[21], vertexData[21]; \
	FieldData data[21]; \
	Arrays backEdgeData_p[21], frontEdgeData_p[21]; \
	Arrays selfConnData_p[21], vertexData_p[21]; \
	FieldData data_p[21]; \
	for(iArray=0;iArray<21;iArray++) \
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
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
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
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
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
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
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
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 4 */ \
	opt = 4; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				phiCDS); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 5 */ \
	opt = 5; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				Xac0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 6 */ \
	opt = 6; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				Xac1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 7 */ \
	opt = 7; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				Xac2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 8 */ \
	opt = 8; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				rcpac); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, rcp); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 9 */ \
	opt = 9; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				denac); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, den); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap; \
	cOpt[opt].fun_host  = swInterpolation_weightedSwap; \
	\
	/* operator 10 */ \
	opt = 10; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 11 */ \
	opt = 11; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 12 */ \
	opt = 12; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 13 */ \
	opt = 13; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 14 */ \
	opt = 14; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 15 */ \
	opt = 15; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 16 */ \
	opt = 16; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				phiUDS); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_directed; \
	cOpt[opt].fun_host  = swInterpolation_directed; \
	\
	/* operator 17 */ \
	opt = 17; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, Xpn0); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_constant_1; \
	cOpt[opt].fun_host  = swInterpolation_constant_1; \
	\
	/* operator 18 */ \
	opt = 18; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, Xpn1); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_constant_1; \
	cOpt[opt].fun_host  = swInterpolation_constant_1; \
	\
	/* operator 19 */ \
	opt = 19; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, Xpn2); \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_constant_1; \
	cOpt[opt].fun_host  = swInterpolation_constant_1; \
	\
	/* operator 20 */ \
	opt = 20; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, phiCDS2);\
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swInterpolation_constant_0_5; \
	cOpt[opt].fun_host  = swInterpolation_constant_0_5; \
	\
}

#define swSeparate_interpolation_test() \
{ \
    int optNum = 21; \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *face_lam, *vis, *visac, *gradPhi0, *dPhidXac0, *gradPhi1; \
	swFloat *dPhidXac1, *gradPhi2, *dPhidXac2; \
	swFloat *mf, *cell_x0, *XU0, *cell_x1, *XU1, *cell_x2, *XU2; \
	swFloat *dPhidXU0, *dPhidXU1, *dPhidXU2, *phi, *phiUDS, *dPhi; \
	swFloat *Xpn0, *Xpn1, *Xpn2, *phiCDS, *Xac0, *Xac1, *Xac2; \
	swFloat *rcp, *rcpac, *den, *denac, *phiCDS2; \
	swFloat *face_lam_m, *vis_m, *visac_m, *gradPhi0_m, *dPhidXac0_m; \
	swFloat *gradPhi1_m, *dPhidXac1_m, *gradPhi2_m, *dPhidXac2_m; \
	swFloat *mf_m, *cell_x0_m, *XU0_m, *cell_x1_m, *XU1_m, *cell_x2_m; \
   	swFloat *XU2_m, *dPhi_m, *phiCDS_m, *Xac0_m, *Xac1_m, *Xac2_m, *rcp_m; \
	swFloat *dPhidXU0_m, *dPhidXU1_m, *dPhidXU2_m, *phi_m, *phiUDS_m; \
	swFloat *Xpn0_m, *Xpn1_m, *Xpn2_m, *rcpac_m, *den_m, *denac_m ; \
	swFloat *phiCDS2_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	constructEmptyArray(paraData); \
	swSeparate_interpolation_data(face_lam, vis, visac, gradPhi0, dPhidXac0, gradPhi1, \
				dPhidXac1, gradPhi2, dPhidXac2, mf, cell_x0, XU0, cell_x1, \
				XU1, cell_x2, XU2, dPhidXU0, dPhidXU1, dPhidXU2, phi, \
				phiUDS, Xpn0, Xpn1, Xpn2, phiCDS, Xac0, Xac1, Xac2, rcp, \
				rcpac, den, denac, phiCDS2, edgeNum, vertexNum); \
	swSeparate_interpolation_data(face_lam_m, vis_m, visac_m, gradPhi0_m, dPhidXac0_m, \
				gradPhi1_m, dPhidXac1_m, gradPhi2_m, dPhidXac2_m, mf_m, \
				cell_x0_m, XU0_m, cell_x1_m, XU1_m, cell_x2_m, XU2_m, \
				dPhidXU0_m, dPhidXU1_m, dPhidXU2_m, phi_m, phiUDS_m, \
				Xpn0_m, Xpn1_m, Xpn2_m, phiCDS_m, Xac0_m, Xac1_m, Xac2_m, \
				rcp_m, rcpac_m, den_m, denac_m, phiCDS2_m, \
				edgeNum, vertexNum); \
	swSeparate_interpolation_master(cOpt_m, face_lam_m, vis_m, visac_m, gradPhi0_m, \
				dPhidXac0_m, gradPhi1_m, dPhidXac1_m, gradPhi2_m, \
				dPhidXac2_m, mf_m, cell_x0_m, XU0_m, cell_x1_m, XU1_m, \
				cell_x2_m, XU2_m, dPhidXU0_m, dPhidXU1_m, dPhidXU2_m, \
				phi_m, phiUDS_m, Xpn0_m, Xpn1_m, Xpn2_m, phiCDS_m, Xac0_m, \
				Xac1_m, Xac2_m, rcp_m, rcpac_m, den_m, denac_m, phiCDS2_m, \
				edgeNum, vertexNum); \
	optNum =21; \
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
	swSeparate_interpolation_host(&mlbIter, face_lam, vis, visac, gradPhi0, dPhidXac0, \
			   	gradPhi1, dPhidXac1, gradPhi2, dPhidXac2, mf, cell_x0, \
				XU0, cell_x1, XU1, cell_x2, XU2,dPhidXU0, dPhidXU1, \
			   	dPhidXU2, phi, phiUDS, Xpn0, Xpn1, Xpn2, phiCDS, Xac0, \
			   	Xac1, Xac2, rcp, rcpac, den, denac, phiCDS2, \
				edgeNum, vertexNum); \
	getTime(time2); \
	printf("Slave Core Time: %f us\n", (time2-time1)*1000000); \
	CG_halt(); \
	\
	checkResult(visac_m, visac, edgeNum); \
	checkResult(rcpac_m, rcpac, edgeNum); \
	checkResult(denac_m, denac, edgeNum); \
	checkResult(phiUDS_m, phiUDS, edgeNum); \
	checkResult(phiCDS_m, phiCDS, edgeNum); \
	checkResult(phiCDS2_m, phiCDS2, edgeNum); \
	checkResult(dPhidXac0_m, dPhidXac0, edgeNum); \
	checkResult(dPhidXac1_m, dPhidXac1, edgeNum); \
	checkResult(dPhidXac2_m, dPhidXac2, edgeNum); \
	checkResult(XU0_m, XU0, edgeNum); \
	checkResult(XU1_m, XU1, edgeNum); \
	checkResult(XU2_m, XU2, edgeNum); \
	checkResult(dPhidXU0_m, dPhidXU0, edgeNum); \
	checkResult(dPhidXU1_m, dPhidXU1, edgeNum); \
	checkResult(dPhidXU2_m, dPhidXU2, edgeNum); \
	checkResult(Xpn0_m, Xpn0, edgeNum); \
	checkResult(Xpn1_m, Xpn1, edgeNum); \
	checkResult(Xpn2_m, Xpn2, edgeNum); \
	checkResult(Xac0_m, Xac0, edgeNum); \
	checkResult(Xac1_m, Xac1, edgeNum); \
	checkResult(Xac2_m, Xac2, edgeNum); \
}

#ifdef __cplusplus
}
#endif

#endif
