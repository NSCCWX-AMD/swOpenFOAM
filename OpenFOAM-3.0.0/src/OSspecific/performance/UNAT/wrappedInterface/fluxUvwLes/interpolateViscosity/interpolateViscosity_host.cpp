#include "interpolateViscosity_host.hpp"
#include "interpolateViscosity.h"

void swInterpolateViscosity_host(MultiLevelBlockIterator *mlbIter,
		   	swFloat *visSx, swFloat *visSy, swFloat *visSz, swFloat *visPNx,
			swFloat *visPNy, swFloat *visPNz, swFloat *rface0, 
			swFloat *rface1, swFloat *face_nx, swFloat *face_ny, 
			swFloat *face_nz, swFloat *face_dx, swFloat *face_dy, 
			swFloat *face_dz, swFloat *face_d, swFloat *viseff, 
			swFloat *facn, swInt edgeNum, swInt vertexNum) 
{ 
	int iArray,opt; 
	coupledOperator *cOpt 
		=(coupledOperator*)malloc(8*sizeof(coupledOperator)); 
	Arrays paraData;
	constructEmptyArray(paraData); 

	Arrays backEdgeData[8], frontEdgeData[8];
	Arrays selfConnData[8], vertexData[8];
	FieldData data[8];
	Arrays backEdgeData_p[8], frontEdgeData_p[8];
	Arrays selfConnData_p[8], vertexData_p[8];
	FieldData data_p[8];
	for(iArray=0;iArray<8;iArray++) 
	{ 
		constructEmptyArray(backEdgeData[iArray]);
		constructEmptyArray(selfConnData[iArray]);
		data[iArray].backEdgeData  = &backEdgeData[iArray];
		data[iArray].frontEdgeData = &frontEdgeData[iArray];
		data[iArray].selfConnData  = &selfConnData[iArray];
		data[iArray].vertexData    = &vertexData[iArray];
		cOpt[iArray].data = &data[iArray];
		
		constructEmptyArray(backEdgeData_p[iArray]);
		constructEmptyArray(selfConnData_p[iArray]);
		data_p[iArray].backEdgeData  = &backEdgeData_p[iArray];
		data_p[iArray].frontEdgeData = &frontEdgeData_p[iArray];
		data_p[iArray].selfConnData  = &selfConnData_p[iArray];
		data_p[iArray].vertexData    = &vertexData_p[iArray];
		cOpt[iArray].data_p = &data_p[iArray];
	} 
	/* operator 0 */ 
	opt = 0;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_nx);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visSx);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, facn);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, COPYIN, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity;
	cOpt[opt].fun_host = swInterpolateViscosity;
	
	/* operator 1 */ 
	opt = 1;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_ny);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visSy);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity;
	cOpt[opt].fun_host = swInterpolateViscosity;
	
	/* operator 2 */ 
	opt = 2;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_nz);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visSz);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity;
	cOpt[opt].fun_host = swInterpolateViscosity;
	
	/* operator 3 */ 
	opt = 3;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_dx);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visPNx);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity;
	cOpt[opt].fun_host = swInterpolateViscosity;
	
	/* operator 4 */ 
	opt = 4;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_dy);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visPNy);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity;
	cOpt[opt].fun_host = swInterpolateViscosity;
	
	/* operator 5 */ 
	opt = 5;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, face_dz);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, visPNz);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity;
	cOpt[opt].fun_host = swInterpolateViscosity;
	
	/* operator 6 */ 
	opt = 6;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface0);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, face_d);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity_1;
	cOpt[opt].fun_host = swInterpolateViscosity_1;
	
	/* operator 7 */ 
	opt = 7;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface1);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, face_d);
	constructSingleArray(vertexData_p[opt], 1, vertexNum, UPDATED, viseff);
	cOpt[opt].fun_slave = slave_swInterpolateViscosity_1;
	cOpt[opt].fun_host = swInterpolateViscosity_1;
	
	mlbIter->edge2VertexIteration(&paraData, cOpt, 8);

	free(cOpt);
	for(iArray=0;iArray<8;iArray++) 
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
}
