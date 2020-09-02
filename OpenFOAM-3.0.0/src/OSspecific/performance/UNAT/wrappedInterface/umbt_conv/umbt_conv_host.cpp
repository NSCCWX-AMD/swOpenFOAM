#include "umbt_conv_host.hpp"
#include "umbt_conv.h"

void umbt_conv_host(Iterator *iter, swFloat *n_vg, swFloat *u, swFloat *cen,
			swFloat *gra, swFloat *rhs, swInt edgeNum,swInt cellNum)
{
	coupledOperator *cOpt=(coupledOperator*)malloc(sizeof(coupledOperator));
	Arrays paraData;
	constructEmptyArray(paraData);

	Arrays frontEdgeData, backEdgeData;
	Arrays selfConnData, vertexData;
	FieldData data;
	constructSingleArray(frontEdgeData, 8, edgeNum, COPYIN, n_vg);
	addSingleArray(frontEdgeData, 3, edgeNum, COPYIN, cen)
	constructEmptyArray(selfConnData);
	constructEmptyArray(backEdgeData);
	constructSingleArray(vertexData, 5, cellNum, COPYOUT, rhs);
	addSingleArray(vertexData, 5, cellNum, COPYIN, u);
	addSingleArray(vertexData, 19, cellNum, COPYIN, gra);
   	data.backEdgeData  = &backEdgeData;
   	data.frontEdgeData = &frontEdgeData;
   	data.selfConnData  = &selfConnData;
   	data.vertexData    = &vertexData;
	cOpt->data = &data;
	cOpt->fun_host = umbt_conv;
	cOpt->fun_slave = slave_umbt_conv;

	//iter->vertex2EdgeIteration(&paraData, cOpt, 1);
	iter->edge2VertexIteration(&paraData, cOpt, 1);

	free(cOpt);
	destroyArray(paraData);
	destroyArray(backEdgeData);
	destroyArray(frontEdgeData);
	destroyArray(selfConnData);
	destroyArray(vertexData);
}


