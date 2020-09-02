#include "calculateFccFlux_host.hpp"
#include "calculateFccFlux.h"

void swCalculateFccFlux_host(MultiLevelBlockIterator *mlbIter, 
		   	swFloat *massFlux, swFloat *fccx, swFloat *dudx, swFloat *dvdx,
		   	swFloat *dwdx, swFloat *Su, swFloat *Sv, swFloat *Sw, 
			swFloat *fccy, swFloat *fccz, swFloat *dudy, swFloat *dvdy, 
			swFloat *dwdy, swFloat *dudz, swFloat *dvdz, swFloat *dwdz, 
			swInt edgeNum, swInt vertexNum) 
{ 
	int iArray,opt; 
	coupledOperator *cOpt 
		=(coupledOperator*)malloc(9*sizeof(coupledOperator)); 
	Arrays paraData; 
	constructEmptyArray(paraData); 

	Arrays backEdgeData[9], frontEdgeData[9]; 
	Arrays selfConnData[9], vertexData[9]; 
	FieldData data[9]; 
	Arrays backEdgeData_p[9], frontEdgeData_p[9]; 
	Arrays selfConnData_p[9], vertexData_p[9]; 
	FieldData data_p[9]; 
	for(iArray=0;iArray<9;iArray++) 
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
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudx);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, COPYIN, fccx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 1 */ 
	opt = 1;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, \
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 2 */ 
	opt = 2;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 3 */ 
	opt = 3;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, COPYIN, fccy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;

	/* operator 4 */ 
	opt = 4;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdy);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 5 */ 
	opt = 5;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 6 */ 
	opt = 6;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, COPYIN, fccz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 7 */ 
	opt = 7;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	/* operator 8 */ 
	opt = 8;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdz);
	addSingleArray(vertexData[opt], 1,vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	addSingleArray(frontEdgeData_p[opt], 1,edgeNum, UPDATED, fccz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateFccFlux;
	cOpt[opt].fun_host  = swCalculateFccFlux;
	
	mlbIter->edge2VertexIteration(&paraData, cOpt, 9);

	free(cOpt);
	destroyArray(paraData);
	for(iArray=0;iArray<9;iArray++)
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

