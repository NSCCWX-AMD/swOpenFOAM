#include "separate_interpolation_host.hpp"
#include "separate_interpolation.h"

void swSeparate_interpolation_host(MultiLevelBlockIterator *mlbIter,
			swFloat* face_lam, swFloat *vis,
			swFloat *visac, swFloat *gradPhi0, swFloat *dPhidXac0, 
			swFloat *gradPhi1, swFloat *dPhidXac1, swFloat *gradPhi2,
			swFloat *dPhidXac2, swFloat *mf, swFloat *cell_x0, 
			swFloat *XU0, swFloat *cell_x1, swFloat *XU1, swFloat* cell_x2,
		   	swFloat *XU2, swFloat *dPhidXU0, swFloat *dPhidXU1, 
			swFloat *dPhidXU2, swFloat *phi, swFloat *phiUDS, swFloat *Xpn0,
		   	swFloat *Xpn1, swFloat *Xpn2, swFloat *phiCDS, swFloat *Xac0, 
			swFloat *Xac1, swFloat *Xac2, swFloat *rcp, swFloat *rcpac, 
			swFloat *den, swFloat *denac, swFloat *phiCDS2, 
			swInt edgeNum, swInt vertexNum) 
{ 
	int iArray,opt; 
	Arrays paraData; 
	constructEmptyArray(paraData); 
	coupledOperator *cOpt 
		=(coupledOperator*)malloc(21*sizeof(coupledOperator)); 
	Arrays backEdgeData[21], frontEdgeData[21];
	Arrays selfConnData[21], vertexData[21];
	FieldData data[21];
	Arrays backEdgeData_p[21], frontEdgeData_p[21];
	Arrays selfConnData_p[21], vertexData_p[21];
	FieldData data_p[21];
	for(iArray=0;iArray<21;iArray++) 
	{ 
		data[iArray].backEdgeData    = &backEdgeData[iArray];
		data[iArray].frontEdgeData   = &frontEdgeData[iArray];
		data[iArray].selfConnData    = &selfConnData[iArray];
		data[iArray].vertexData      = &vertexData[iArray];
		cOpt[iArray].data = &data[iArray];
		data_p[iArray].backEdgeData  = &backEdgeData_p[iArray];
		data_p[iArray].frontEdgeData = &frontEdgeData_p[iArray];
		data_p[iArray].selfConnData  = &selfConnData_p[iArray];
		data_p[iArray].vertexData    = &vertexData_p[iArray];
		cOpt[iArray].data_p = &data_p[iArray];
	} 
	/* operator 0 */ 
	opt = 0;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visac);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, vis);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 1 */ 
	opt = 1;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXac0);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 2 */ 
	opt = 2;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXac1);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 3 */ 
	opt = 3;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXac2);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 4 */ 
	opt = 4;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				phiCDS);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 5 */ 
	opt = 5;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				Xac0);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 6 */ 
	opt = 6;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				Xac1);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 7 */ 
	opt = 7;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				Xac2);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 8 */ 
	opt = 8;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				rcpac);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, rcp);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 9 */ 
	opt = 9;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				denac);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, den);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				face_lam);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_weightedSwap;
	cOpt[opt].fun_host  = swInterpolation_weightedSwap;
	
	/* operator 10 */ 
	opt = 10;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU0);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 11 */ 
	opt = 11;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU1);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 12 */ 
	opt = 12;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				XU2);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 13 */ 
	opt = 13;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU0);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi0);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 14 */ 
	opt = 14;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU1);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi1);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 15 */ 
	opt = 15;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				dPhidXU2);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, gradPhi2);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 16 */ 
	opt = 16;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, \
				phiUDS);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, mf);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_directed;
	cOpt[opt].fun_host  = swInterpolation_directed;
	
	/* operator 17 */ 
	opt = 17;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, Xpn0);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x0);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructEmptyArray(frontEdgeData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_constant_1;
	cOpt[opt].fun_host  = swInterpolation_constant_1;
	
	/* operator 18 */ 
	opt = 18;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, Xpn1);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x1);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructEmptyArray(frontEdgeData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_constant_1;
	cOpt[opt].fun_host  = swInterpolation_constant_1;
	
	/* operator 19 */ 
	opt = 19;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, Xpn2);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cell_x2);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructEmptyArray(frontEdgeData_p[opt]);
	constructEmptyArray(vertexData_p[opt]);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_constant_1;
	cOpt[opt].fun_host  = swInterpolation_constant_1;
	
	/* operator 20 */ 
	opt = 20;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, phiCDS2);\
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, phi);
	constructEmptyArray(backEdgeData[opt]);
	constructEmptyArray(selfConnData[opt]);
	constructEmptyArray(vertexData_p[opt]);
	constructEmptyArray(backEdgeData_p[opt]);
	constructEmptyArray(frontEdgeData_p[opt]);
	constructEmptyArray(selfConnData_p[opt]);
	cOpt[opt].fun_slave = slave_swInterpolation_constant_0_5;
	cOpt[opt].fun_host  = swInterpolation_constant_0_5;
	
	mlbIter->edge2VertexIteration(&paraData, cOpt, 21);

	for(iArray=0;iArray<21;iArray++) 
	{ 
		destroyArray(frontEdgeData[iArray]);
		destroyArray(backEdgeData[iArray]);
		destroyArray(selfConnData[iArray]);
		destroyArray(vertexData[iArray]);

		destroyArray(frontEdgeData_p[iArray]);
		destroyArray(backEdgeData_p[iArray]);
		destroyArray(selfConnData_p[iArray]);
		destroyArray(vertexData_p[iArray]);
	} 
	free(cOpt);
}

