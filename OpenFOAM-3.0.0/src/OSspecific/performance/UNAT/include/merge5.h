#ifndef MERGE5_H
#define MERGE5_H
#include "swMacro.h"
#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

define_e2v_FunPtr(swMerge5_mf_1);
define_e2v_slaveFunPtr(swMerge5_mf_1);
define_e2v_FunPtr(swMerge5_mf_2);
define_e2v_slaveFunPtr(swMerge5_mf_2);
define_e2v_FunPtr(swMerge5_mf_3);
define_e2v_slaveFunPtr(swMerge5_mf_3);
define_e2v_FunPtr(swMerge5_lam_1);
define_e2v_slaveFunPtr(swMerge5_lam_1);
define_e2v_FunPtr(swMerge5_lam_2);
define_e2v_slaveFunPtr(swMerge5_lam_2);
define_e2v_FunPtr(swMerge5_lam_3);
define_e2v_slaveFunPtr(swMerge5_lam_3);
define_e2v_FunPtr(swMerge5_lam_4);
define_e2v_slaveFunPtr(swMerge5_lam_4);
define_e2v_FunPtr(swMerge5_su_1);
define_e2v_slaveFunPtr(swMerge5_su_1);
define_e2v_FunPtr(swMerge5_su_2);
define_e2v_slaveFunPtr(swMerge5_su_2);
define_e2v_FunPtr(swMerge5_su_3);
define_e2v_slaveFunPtr(swMerge5_su_3);
define_e2v_FunPtr(swMerge5_su_4);
define_e2v_slaveFunPtr(swMerge5_su_4);

#define swMerge5_data(mf, cell_x0, face_x0, XacU0, cell_x1, face_x1, \
			XacU1, cell_x2, face_x2, XacU2, graphPhi0, phiDelta, graphPhi1,\
		    graphPhi2, phi, su, \
			face_lam, vis, face_D, visFace, visac, dPhi, dPhidFn, face_n0,\
			face_n1, face_n2, rface1, rface2, fdi, fde,edgeNum, vertexNum) \
{ \
	int i,j; \
	mf       = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x0  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x1  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_x2  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XacU0    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XacU1    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	XacU2    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	phiDelta = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_lam = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_D   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visFace  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	visac    = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhi     = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	dPhidFn  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n0  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n1  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	face_n2  = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface1   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	rface2   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fdi      = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	fde   = (swFloat*)malloc(edgeNum*sizeof(swFloat)); \
	swFloat* edge[21] = {mf, face_x0, face_x1, face_x2, XacU0, XacU1, \
		XacU2, phiDelta, face_lam, face_D, visFace, visac, dPhi, dPhidFn, \
		face_n0, face_n1, face_n2, rface1, rface2, fdi, fde} ;\
	for(i=0;i<21;i++) \
	{ \
		for(j=0;j<edgeNum;j++) \
		{ \
			edge[i][j] = (swFloat)(j+i+1)/(j+i+2); \
		} \
	} \
	cell_x0  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cell_x1  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	cell_x2  = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	gradPhi0 = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	gradPhi1 = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	gradPhi2 = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	phi      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	su       = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	vis      = (swFloat*)malloc(vertexNum*sizeof(swFloat)); \
	swFloat *cell[9] = {cell_x0, cell_x1, cell_x2, gradPhi0, gradPhi1, \
		gradPhi2, phi, su, vis}; \
	for(i=0;i<9;i++) \
	{\
		for(j=0;j<vertexNum;j++) \
		{ \
			cell[i][j] = (swFloat)(2*j+i+1)/(2*j+i+2); \
		} \
	} \
} 

#define swMerge5_host(cOpt, mf, cell_x0, face_x0, XacU0, cell_x1, \
			face_x1, XacU1, cell_x2, face_x2, XacU2, graphPhi0, phiDelta, \
			graphPhi1, graphPhi2, phi, su, \
			face_lam, vis, face_D, visFace, visac, dPhi, dPhidFn, face_n0, \
			face_n1, face_n2, rface1, rface2, fdi, fde,edgeNum, vertexNum) \
{ \
	int iArray,opt,iOpt; \
	Arrays backEdgeData[21], frontEdgeData[21]; \
	Arrays selfConnData[21], vertexData[21]; \
	Arrays backEdgeData_p[21], frontEdgeData_p[21]; \
	Arrays selfConnData_p[21], vertexData_p[21]; \
	FieldData data[21], data_p[21]; \
	for(iOpt=0;iOpt<21;iOpt++) \
	{ \
		data[iOpt].backEdgeData  = &backEdgeData[iOpt]; \
		data[iOpt].frontEdgeData = &frontEdgeData[iOpt]; \
		data[iOpt].selfConnData  = &selfConnData[iOpt]; \
		data[iOpt].vertexData    = &vertexData[iOpt]; \
		data_p[iOpt].backEdgeData  = &backEdgeData_p[iOpt]; \
		data_p[iOpt].frontEdgeData = &frontEdgeData_p[iOpt]; \
		data_p[iOpt].selfConnData  = &selfConnData_p[iOpt]; \
		data_p[iOpt].vertexData    = &vertexData_p[iOpt]; \
		cOpt[iOpt].data   = &data[iOpt]; \
		cOpt[iOpt].data_p = &data_p[iOpt]; \
	} \
	\
	/* operator 0 */ \
	opt = 0; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, XacU0); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  face_x0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_1; \
	cOpt[opt].fun_host  = swMerge5_mf_1; \
	\
	/* operator 1 */ \
	opt = 1; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, XacU1); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  face_x1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_1; \
	cOpt[opt].fun_host  = swMerge5_mf_1; \
	\
	/* operator 2 */ \
	opt = 2; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	constructSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, XacU2); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  face_x2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_1; \
	cOpt[opt].fun_host  = swMerge5_mf_1; \
	\
	/* operator 3 */ \
	opt = 3; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0); \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  XacU0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	addSingleArray(frontEdgeData_p[opt],1,edgeNum, COPYIN, phiDelta); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_2; \
	cOpt[opt].fun_host  = swMerge5_mf_2; \
	\
	/* operator 4 */ \
	opt = 4; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1); \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  XacU1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	addSingleArray(frontEdgeData_p[opt],1,edgeNum, UPDATED, phiDelta); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_2; \
	cOpt[opt].fun_host  = swMerge5_mf_2; \
	\
	/* operator 5 */ \
	opt = 5; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2); \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  XacU2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	addSingleArray(frontEdgeData_p[opt],1,edgeNum, UPDATED, phiDelta); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_2; \
	cOpt[opt].fun_host  = swMerge5_mf_2; \
	\
	/* operator 6 */ \
	opt = 6; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT,  su); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	addSingleArray(frontEdgeData_p[opt],1,edgeNum, COPYOUT, phiDelta); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_mf_3; \
	cOpt[opt].fun_host  = swMerge5_mf_3; \
	\
	/* operator 7 */ \
	opt = 7; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN,  face_D); \
	addSingleArray(frontEdgeData[opt],1,edgeNum, COPYOUT, visFace); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				face_lam); \
	constructSingleArray(vertexData_p[opt], 1, vertexNum, COPYIN, vis); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_1; \
	cOpt[opt].fun_host  = swMerge5_lam_1; \
	\
	/* operator 8 */ \
	opt = 8; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT,  visac); \
	constructEmptyArray(vertexData[opt]); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, vis); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_2; \
	cOpt[opt].fun_host  = swMerge5_lam_2; \
	\
	/* operator 9 */ \
	opt = 9; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYINOUT,  dPhi); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_3; \
	cOpt[opt].fun_host  = swMerge5_lam_3; \
	\
	/* operator 10 */ \
	opt = 10; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED,  dPhi); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_3; \
	cOpt[opt].fun_host  = swMerge5_lam_3; \
	\
	/* operator 11 */ \
	opt = 11; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2); \
	addSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2); \
	constructEmptyArray(frontEdgeData[opt]); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYOUT,  dPhi); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_3; \
	cOpt[opt].fun_host  = swMerge5_lam_3; \
	\
	/* operator 12 */ \
	opt = 12; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0); \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_n0); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN,  dPhidFn); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_4; \
	cOpt[opt].fun_host  = swMerge5_lam_4; \
	\
	/* operator 13 */ \
	opt = 13; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1); \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_n1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED,  dPhidFn); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_4; \
	cOpt[opt].fun_host  = swMerge5_lam_4; \
	\
	/* operator 14 */ \
	opt = 14; \
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2); \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_n2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYOUT,  dPhidFn); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_lam_4; \
	cOpt[opt].fun_host  = swMerge5_lam_4; \
	\
	/* operator 15 */ \
	opt = 15; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, visac); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidFn); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fde); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_su_4; \
	cOpt[opt].fun_host  = swMerge5_su_4; \
	\
	/* operator 16 */ \
	opt = 16; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, visFace); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhi); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fdi); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_su_4; \
	cOpt[opt].fun_host  = swMerge5_su_4; \
	\
	/* operator 17 */ \
	opt = 17; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, visac); \
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhidFn); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(vertexData_p[opt], 1, vertexNum, COPYOUT, su); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(frontEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_su_2; \
	cOpt[opt].fun_host  = swMerge5_su_2; \
	\
	/* operator 18 */ \
	opt = 18; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, dPhi); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(vertexData_p[opt], 1, vertexNum, COPYOUT, su); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				visFace); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_su_1; \
	cOpt[opt].fun_host  = swMerge5_su_1; \
	\
	/* operator 19 */ \
	opt = 19; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, rface1); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				visFace); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_su_3; \
	cOpt[opt].fun_host  = swMerge5_su_3; \
	\
	/* operator 20 */ \
	opt = 20; \
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, rface2); \
	constructEmptyArray(backEdgeData[opt]); \
	constructEmptyArray(selfConnData[opt]); \
	constructEmptyArray(vertexData[opt]); \
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				visFace); \
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf); \
	constructEmptyArray(backEdgeData_p[opt]); \
	constructEmptyArray(selfConnData_p[opt]); \
	constructEmptyArray(vertexData_p[opt]); \
	cOpt[opt].fun_slave = slave_swMerge5_su_3; \
	cOpt[opt].fun_host  = swMerge5_su_3; \
	\
}

#define swMerge5_test() \
{ \
    int optNum = 21; \
	coupledOperator *cOpt \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	coupledOperator *cOpt_m \
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator)); \
	swFloat *mf, *cell_x0, *cell_x1, *cell_x2, *face_x0, *face_x1; \
	swFloat *face_x2, *XacU0, *XacU1, *XacU2, *gradPhi0, *phiDelta; \
	swFloat *gradPhi1, *gradPhi2, *phi, *su, *fde, *fdi; \
	swFloat *face_lam, *vis, *face_D, *visFace, *visac, *dPhi; \
	swFloat *dPhidFn, *face_n0, *face_n1, *face_n2, *rface1, *rface2; \
	swFloat *mf_m, *cell_x0_m, *cell_x1_m, *cell_x2_m, *face_x0_m; \
	swFloat *face_x1_m, *face_x2_m, *XacU0_m, *XacU1_m, *XacU2_m; \
	swFloat *gradPhi0_m, *phiDelta_m, *gradPhi1_m, *gradPhi2_m; \
	swFloat *phi_m, *su_m, *rface2_m, *fde_m, *fdi_m; \
	swFloat *face_lam_m, *vis_m, *face_D_m, *visFace_m, *visac_m, *dPhi_m; \
	swFloat *dPhidFn_m, *face_n0_m, *face_n1_m, *face_n2_m, *rface1_m; \
	Arrays backEdgeData_c, frontEdgeData_c, selfConnData_c, vertexData_c; \
	Arrays backEdgeData_m, frontEdgeData_m, selfConnData_m, vertexData_m; \
	Arrays paraData; \
	swFloat paras[4] = {1.111111, 2.2222222, 3.333333, 4.4444444}; \
	constructSingleArray(paraData, 1, 4, COPYIN, &paras[0]); \
	swMerge5_data(mf, cell_x0, face_x0, XacU0, cell_x1, face_x1, \
				XacU1, cell_x2, face_x2, XacU2, gradPhi0, phiDelta, \
				gradPhi1, gradPhi2, phi, su, face_lam, vis, face_D, \
				visFace, visac, dPhi, dPhidFn, face_n0, face_n1, face_n2, \
				rface1, rface2, fdi, fde, edgeNum, vertexNum); \
	swMerge5_host(cOpt, mf, cell_x0, face_x0, XacU0, cell_x1, face_x1, \
				XacU1, cell_x2, face_x2, XacU2, gradPhi0, phiDelta, \
				gradPhi1, gradPhi2, phi, su, face_lam, vis, face_D, \
				visFace, visac, dPhi, dPhidFn, face_n0, face_n1, face_n2, \
				rface1, rface2, fdi, fde, edgeNum, vertexNum); \
	swMerge5_data(mf_m, cell_x0_m, face_x0_m, XacU0_m, cell_x1_m, \
				face_x1_m, XacU1_m, cell_x2_m, face_x2_m, XacU2_m, \
				gradPhi0_m, phiDelta_m, gradPhi1_m, gradPhi2_m, phi_m, \
				su_m, face_lam_m, vis_m, face_D_m, visFace_m, visac_m, \
				dPhi_m, dPhidFn_m, face_n0_m, face_n1_m, face_n2_m, \
				rface1_m, rface2_m, fdi_m, fde_m, edgeNum, vertexNum); \
	swMerge5_host(cOpt_m, mf_m, cell_x0_m, face_x0_m, XacU0_m, \
				cell_x1_m, face_x1_m, XacU1_m, cell_x2_m, face_x2_m, \
				XacU2_m, gradPhi0_m, phiDelta_m, gradPhi1_m, gradPhi2_m, \
				phi_m, su_m, face_lam_m, vis_m, face_D_m, visFace_m, \
				visac_m, dPhi_m, dPhidFn_m, face_n0_m, face_n1_m, \
				face_n2_m, rface1_m, rface2_m, fdi_m, fde_m, \
				edgeNum, vertexNum); \
	\
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
	printf("test\n");\
	optNum = 17;\
	getTime(time1); \
	for(int iOpt=0;iOpt<optNum;iOpt++) \
	{ \
		backEdgeData_m = *(cOpt_m[iOpt].data_p->backEdgeData); \
		frontEdgeData_m = *(cOpt_m[iOpt].data_p->frontEdgeData); \
		selfConnData_m = *(cOpt_m[iOpt].data_p->selfConnData); \
		vertexData_m = *(cOpt_m[iOpt].data_p->vertexData); \
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
	checkResult(XacU0_m, XacU0, edgeNum*dims); \
	checkResult(XacU1_m, XacU1, edgeNum*dims); \
	checkResult(XacU2_m, XacU2, edgeNum*dims); \
	checkResult(phiDelta_m, phiDelta, edgeNum*dims); \
	checkResult(su_m, su, vertexNum*dims); \
	checkResult(visFace_m, visFace, edgeNum*dims); \
	checkResult(visac_m, visac, edgeNum*dims); \
	checkResult(dPhi_m, dPhi, edgeNum*dims); \
	checkResult(dPhidFn_m, dPhidFn, edgeNum*dims); \
	checkResult(rface1_m, rface1, edgeNum*dims); \
	checkResult(rface2_m, rface2, edgeNum*dims); \
}

#ifdef __cplusplus
}
#endif

#endif
