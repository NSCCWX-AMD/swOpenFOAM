#include "separate_integrate_host.hpp"
#include "separate_integrate.h"

void swSeparateIntegrate_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *fu, swFloat *su, swInt edgeNum, swInt vertexNum) 
{ 
	int iArray,opt; 
	coupledOperator *cOpt
		= (coupledOperator*)malloc(sizeof(coupledOperator));
	Arrays paraData;
	constructEmptyArray(paraData);

	Arrays backEdgeData[1], frontEdgeData[1]; 
	Arrays selfConnData[1], vertexData[1]; 
	FieldData data[1]; 
	Arrays backEdgeData_p[1], frontEdgeData_p[1]; 
	Arrays selfConnData_p[1], vertexData_p[1]; 
	FieldData data_p[1]; 
	for(iArray=0;iArray<1;iArray++) 
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
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, fu); 
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, su);
	constructEmptyArray(backEdgeData[opt]); 
	constructEmptyArray(selfConnData[opt]); 
	constructEmptyArray(frontEdgeData_p[opt]); 
	constructEmptyArray(backEdgeData_p[opt]); 
	constructEmptyArray(selfConnData_p[opt]); 
	constructEmptyArray(vertexData_p[opt]); 
	cOpt[opt].fun_slave = slave_swSeparateIntegrate_su; 
	cOpt[opt].fun_host  = swSeparateIntegrate_su; 

	mlbIter->edge2VertexIteration(&paraData, cOpt, 1);

	destroyArray(frontEdgeData[0]);
	destroyArray(backEdgeData[0]);
	destroyArray(selfConnData[0]);
	destroyArray(vertexData[0]);
	destroyArray(frontEdgeData_p[0]);
	destroyArray(backEdgeData_p[0]);
	destroyArray(selfConnData_p[0]);
	destroyArray(vertexData_p[0]);
	destroyArray(paraData);
	free(cOpt);
}
