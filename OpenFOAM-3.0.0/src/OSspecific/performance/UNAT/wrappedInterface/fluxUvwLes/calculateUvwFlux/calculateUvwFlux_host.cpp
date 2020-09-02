#include "calculateUvwFlux_host.hpp"
#include "calculateUvwFlux.h"

void swCalculateUvwFlux_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *massFlux, swFloat *u, swFloat *v, swFloat *w, 
			swFloat *Su, swFloat *Sv, swFloat *Sw, swFloat gamblend, 
			swInt edgeNum, swInt vertexNum) 
{ 
	int iArray,opt;
	coupledOperator *cOpt 
		=(coupledOperator*)malloc(3*sizeof(coupledOperator));
	Arrays paraData;
	swFloat paras[1] = {gamblend};
	constructSingleArray(paraData, 1, 1, COPYIN, &paras[0]);

	Arrays backEdgeData[3], frontEdgeData[3];
	Arrays selfConnData[3], vertexData[3];
	FieldData data[3];
	Arrays backEdgeData_p[3], frontEdgeData_p[3];
	Arrays selfConnData_p[3], vertexData_p[3];
	FieldData data_p[3];
	for(iArray=0;iArray<3;iArray++) 
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
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, u);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_host  = swCalculateUvwFlux;
	cOpt[opt].fun_slave = slave_swCalculateUvwFlux;
	
	/* operator 1 */ 
	opt = 1;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, v);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_host  = swCalculateUvwFlux;
	cOpt[opt].fun_slave = slave_swCalculateUvwFlux;
	
	/* operator 2 */ 
	opt = 2;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, w);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_host  = swCalculateUvwFlux;
	cOpt[opt].fun_slave = slave_swCalculateUvwFlux;
	
	mlbIter->edge2VertexIteration(&paraData, cOpt, 3);

	free(cOpt);
	destroyArray(paraData);
	for(iArray=0;iArray<3;iArray++)
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

