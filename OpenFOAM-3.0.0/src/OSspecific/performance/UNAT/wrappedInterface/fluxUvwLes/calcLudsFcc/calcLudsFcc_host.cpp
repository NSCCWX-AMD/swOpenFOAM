#include "calcLudsFcc_host.hpp"
#include "calcLudsFcc.h"

void swCalcLudsFcc_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *massFlux, swFloat *facex_x, swFloat *facex_y, 
			swFloat *facex_z, swFloat *cellx_x, swFloat *cellx_y, 
			swFloat *cellx_z, swFloat *fccx, swFloat *fccy, swFloat *fccz,
			swFloat *rface0, swFloat *rface1, swFloat *facn, swFloat CS, 
			swFloat gamblend, swInt edgeNum, swInt vertexNum) 
{ 
	int iArray,opt,optNum; 
	if(CS == 0 || CS == 1)
	{
		optNum = 1;
	} else if(CS == 3)
	{
		optNum = 4;
	} else
	{
		//pr_error("!!!");
	}
	coupledOperator *cOpt  
		=(coupledOperator*)malloc(optNum*sizeof(coupledOperator));
	swFloat paras[2] = {gamblend, CS}; 
	Arrays paraData;
	constructSingleArray(paraData, 1, 2, COPYIN, &paras[0]); 

	Arrays backEdgeData[optNum], frontEdgeData[optNum]; 
	Arrays selfConnData[optNum], vertexData[optNum];
	FieldData data[optNum];
	Arrays backEdgeData_p[optNum], frontEdgeData_p[optNum];
	Arrays selfConnData_p[optNum], vertexData_p[optNum];
	FieldData data_p[optNum];
	for(iArray=0;iArray<optNum;iArray++) 
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
	if(CS == 3)
	{
	/* operator 0 */ 
	opt = 0;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fccx);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, facex_x);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cellx_x);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalcLudsFcc;
	cOpt[opt].fun_host  = swCalcLudsFcc;
	
	/* operator 1 */ 
	opt = 1;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fccy);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, facex_y);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cellx_y);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalcLudsFcc;
	cOpt[opt].fun_host  = swCalcLudsFcc;
	
	/* operator 2 */ 
	opt = 2;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum, COPYOUT, fccz);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYIN, facex_z);
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, cellx_z);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalcLudsFcc;
	cOpt[opt].fun_host  = swCalcLudsFcc;
	
	/* operator 3 */ 
	opt = 3;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface0);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYINOUT, rface1);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalcLudsFcc_rface;
	cOpt[opt].fun_host  = swCalcLudsFcc_rface;
	} else if(CS == 1)
	{
	/* operator 0 */ 
	opt = 0;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface0);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYINOUT, rface1);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, 
				facn);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalcLudsFcc_rface_CDS;
	cOpt[opt].fun_host  = swCalcLudsFcc_rface_CDS;
	} else if(CS == 0)
	{
	/* operator 0 */ 
	opt = 0;
	constructSingleArray(frontEdgeData[opt], 1, edgeNum,COPYINOUT,rface0);
	addSingleArray(frontEdgeData[opt], 1, edgeNum, COPYINOUT, rface1);
	constructEmptyArray(vertexData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, 
				massFlux);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalcLudsFcc_rface;
	cOpt[opt].fun_host  = swCalcLudsFcc_rface;
	}
	
	mlbIter->edge2VertexIteration(&paraData, cOpt, optNum);

	free(cOpt);
	destroyArray(paraData);
	for(iArray=0;iArray<optNum;iArray++)
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

