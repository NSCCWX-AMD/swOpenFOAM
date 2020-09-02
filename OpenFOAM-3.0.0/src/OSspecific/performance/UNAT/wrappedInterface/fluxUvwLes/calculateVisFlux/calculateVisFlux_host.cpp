#include "calculateVisFlux_host.hpp"
#include "calculateVisFlux.h"

void swCalculateVisFlux_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *facn, swFloat *visSy, swFloat *dudy, swFloat *dvdx, 
			swFloat *dvdy, swFloat *dvdz, swFloat *dwdy, swFloat *Su, 
			swFloat *Sv, swFloat *Sw, swFloat *dudx, swFloat *dudz, 
			swFloat *dwdx, swFloat *dwdz, swFloat *visSx, swFloat *visSz,
			swFloat *visPNx, swFloat *visPNy, swFloat *visPNz, 
			swInt edgeNum, swInt vertexNum) 
{
	int iArray,opt; 
	coupledOperator *cOpt 
		=(coupledOperator*)malloc(24*sizeof(coupledOperator)); 
	Arrays paraData; 
	constructEmptyArray(paraData); 

	Arrays backEdgeData[24], frontEdgeData[24];
	Arrays selfConnData[24], vertexData[24];
	FieldData data[24];
	Arrays backEdgeData_p[24], frontEdgeData_p[24];
	Arrays selfConnData_p[24], vertexData_p[24];
	FieldData data_p[24];
	for(iArray=0;iArray<24;iArray++) 
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
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visSy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 1 */ 
	opt = 1;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 2 */ 
	opt = 2;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux_multiply2;
	cOpt[opt].fun_host  = swCalculateVisFlux_multiply2;
	
	/* operator 3 */ 
	opt = 3;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum,UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum,UPDATED, visSy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 4 */ 
	opt = 4;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 5 */ 
	opt = 5;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visSx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux_multiply2;
	cOpt[opt].fun_host  = swCalculateVisFlux_multiply2;
	
	/* operator 6 */ 
	opt = 6;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 7 */ 
	opt = 7;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 8 */ 
	opt = 8;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 9 */ 
	opt = 9;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 10 */ 
	opt = 10;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visSz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 11 */ 
	opt = 11;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 12 */ 
	opt = 12;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 13 */ 
	opt = 13;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 14 */ 
	opt = 14;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visSz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux_multiply2;
	cOpt[opt].fun_host  = swCalculateVisFlux_multiply2;
	
	/* operator 15 */ 
	opt = 15;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visPNx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 16 */ 
	opt = 16;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 17 */ 
	opt = 17;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdx);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNx);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 18 */ 
	opt = 18;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visPNy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 19 */ 
	opt = 19;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 20 */ 
	opt = 20;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdy);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNy);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 21 */ 
	opt = 21;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dudz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Su);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, COPYIN, visPNz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 22 */ 
	opt = 22;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dvdz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sv);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	/* operator 23 */ 
	opt = 23;
	constructSingleArray(vertexData[opt], 1, vertexNum, COPYIN, dwdz);
	addSingleArray(vertexData[opt], 1, vertexNum, COPYOUT, Sw);
	constructEmptyArray(frontEdgeData[opt]);
	constructSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, facn);
	addSingleArray(frontEdgeData_p[opt], 1, edgeNum, UPDATED, visPNz);
	constructEmptyArray(vertexData_p[opt]);
	cOpt[opt].fun_slave = slave_swCalculateVisFlux;
	cOpt[opt].fun_host  = swCalculateVisFlux;
	
	mlbIter->edge2VertexIteration(&paraData, cOpt, 24);

	free(cOpt);
	destroyArray(paraData);
	for(iArray=0;iArray<24;iArray++)
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

