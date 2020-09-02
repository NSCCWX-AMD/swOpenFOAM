#include "spMV_host.hpp"
#include "spMV.h"

void swSpMV_v2e_host(Iterator *iter, swFloat *edge, 
			swFloat *bPtr, swFloat *xPtr, swFloat *diagPtr, 
			swInt dims, swInt edgeNum,swInt cellNum)
{
	int optNum = 2;
	coupledOperator *cOpt
		= (coupledOperator*)malloc(optNum*sizeof(coupledOperator));
	Arrays paraData;
	constructEmptyArray(paraData);

	Arrays frontEdgeData[optNum], backEdgeData[optNum];
	Arrays selfConnData[optNum], vertexData[optNum];
	FieldData data[optNum];
	constructSingleArray(frontEdgeData[0], dims, edgeNum, COPYIN, edge);
	constructSingleArray(selfConnData[0], dims, cellNum, COPYIN, diagPtr);
	constructSingleArray(vertexData[0], dims, cellNum, COPYIN, xPtr);
	addSingleArray(vertexData[0], dims, cellNum, COPYOUT, bPtr);
	constructEmptyArray(backEdgeData[0]);
   	data[0].backEdgeData  = &backEdgeData[0];
   	data[0].frontEdgeData = &frontEdgeData[0];
   	data[0].selfConnData  = &selfConnData[0];
   	data[0].vertexData    = &vertexData[0];
	cOpt[0].data = &data[0];
	cOpt[0].fun_host = swSpMV;
	cOpt[0].fun_slave = slave_swSpMV;

	constructSingleArray(frontEdgeData[1], dims, edgeNum, UPDATED, edge);
	constructSingleArray(selfConnData[1], dims, cellNum, COPYIN, diagPtr);
	constructSingleArray(vertexData[1], dims, cellNum, UPDATED, xPtr);
	addSingleArray(vertexData[1], dims, cellNum, COPYOUT, bPtr);
	constructEmptyArray(backEdgeData[1]);
   	data[1].backEdgeData  = &backEdgeData[1];
   	data[1].frontEdgeData = &frontEdgeData[1];
   	data[1].selfConnData  = &selfConnData[1];
   	data[1].vertexData    = &vertexData[1];
	cOpt[1].data = &data[1];
	cOpt[1].fun_host = swSpMV;
	cOpt[1].fun_slave = slave_swSpMV;

	optNum = 1;
	iter->vertex2EdgeIteration(&paraData, cOpt, optNum);

	free(cOpt);
	destroyArray(paraData);
	for(int iOpt=0;iOpt<optNum;iOpt++)
	{
		destroyArray(backEdgeData[iOpt]);
		destroyArray(frontEdgeData[iOpt]);
		destroyArray(selfConnData[iOpt]);
		destroyArray(vertexData[iOpt]);
	}
}

void swSpMV_e2v_host(Iterator *iter, swFloat *lowerPtr, swFloat *upperPtr, 
			swFloat *bPtr, swFloat *xPtr, swFloat *diagPtr,
			swInt dims, swInt edgeNum,swInt cellNum)
{
	int optNum = 1;
	int iOpt;
	coupledOperator *cOpt
		= (coupledOperator*)malloc(optNum*sizeof(coupledOperator));
	Arrays paraData;
	constructEmptyArray(paraData);

	Arrays frontEdgeData[optNum], backEdgeData[optNum];
	Arrays selfConnData[optNum], vertexData[optNum];
	Arrays frontEdgeData_p[optNum], backEdgeData_p[optNum];
	Arrays selfConnData_p[optNum], vertexData_p[optNum];

	FieldData data[optNum], data_p[optNum];

	for(iOpt=0;iOpt<optNum;iOpt++)
	{
		data[iOpt].backEdgeData  = &backEdgeData[iOpt];
		data[iOpt].frontEdgeData = &frontEdgeData[iOpt];
		data[iOpt].selfConnData  = &selfConnData[iOpt];
		data[iOpt].vertexData    = &vertexData[iOpt];
		cOpt[iOpt].data = &data[iOpt];

		constructEmptyArray(frontEdgeData_p[iOpt]);
		constructEmptyArray(backEdgeData_p[iOpt]);
		constructEmptyArray(selfConnData_p[iOpt]);
		constructEmptyArray(vertexData_p[iOpt]);
		data_p[iOpt].backEdgeData  = &backEdgeData_p[iOpt];
		data_p[iOpt].frontEdgeData = &frontEdgeData_p[iOpt];
		data_p[iOpt].selfConnData  = &selfConnData_p[iOpt];
		data_p[iOpt].vertexData    = &vertexData_p[iOpt];
		cOpt[iOpt].data_p = &data_p[iOpt];
	}

	iOpt = 0;
	constructSingleArray(frontEdgeData[iOpt], dims, edgeNum, COPYIN,
				upperPtr);
	constructSingleArray(selfConnData[iOpt], dims, cellNum, COPYIN, 
				diagPtr);
	constructSingleArray(vertexData[iOpt], dims, cellNum, COPYIN, xPtr);
	addSingleArray(vertexData[iOpt], dims, cellNum, COPYOUT, bPtr);
	constructSingleArray(backEdgeData[iOpt], dims, edgeNum, COPYIN,
				lowerPtr);
	cOpt[iOpt].fun_host = swSpMV;
	cOpt[iOpt].fun_slave = slave_swSpMV;

	iter->edge2VertexIteration(&paraData, cOpt, optNum);

	free(cOpt);
	destroyArray(paraData);
	for(int iOpt=0;iOpt<optNum;iOpt++)
	{
		destroyArray(backEdgeData[iOpt]);
		destroyArray(frontEdgeData[iOpt]);
		destroyArray(selfConnData[iOpt]);
		destroyArray(vertexData[iOpt]);
	}
}


