#include "separate_vector_host.hpp"
#include "separate_vector.h"

void swSeparate_vector_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *mf, swFloat *phiUDS, swFloat *face_x0, 
			swFloat *dPhidXac0, swFloat *face_x1, swFloat *dPhidXac1,
			swFloat *face_x2, swFloat *dPhidXac2, swFloat *XU0, 
			swFloat *XU1, swFloat *XU2, swFloat *visac, swFloat *face_D,
			swFloat *Xpn0, swFloat *Xpn1, swFloat *Xpn2, swFloat *face_n0, 
			swFloat *face_n1, swFloat *face_n2, swFloat *dPhidXU0, 
			swFloat *dPhidXU1, swFloat *dPhidXU2, swFloat *fu, 
			swFloat *rface_1, swFloat *rface_2, swFloat *rcpac, 
			swFloat *denac, swFloat *Xac0, swFloat *Xac1, swFloat *Xac2, 
			swFloat *phiCDS, swFloat *phiCDS2, swFloat csBlend,
			swFloat vis_lam, swFloat vis_sigma, swFloat vis_PrScNr, 
			swFloat vis_lambda, swFloat vis_diffuse, swFloat vis_type,
			swFloat cvScheme, swFloat csFaceCorrect,
			swInt edgeNum, swInt vertexNum) 
{ 
	int iArray; 
	coupledArrayOperator *cOpt 
		=(coupledArrayOperator*)malloc(sizeof(coupledArrayOperator));
	Arrays paraData;
	swFloat paras[9] = {csBlend, vis_lam, vis_sigma, vis_PrScNr, 
		vis_lambda, vis_diffuse, vis_type, cvScheme, 
		csFaceCorrect};
	constructSingleArray(paraData,1,9,COPYIN,&paras[0]); 

	Arrays backEdgeData, frontEdgeData; 
	Arrays selfConnData, vertexData; 
	FieldData data; 
	constructSingleArray(frontEdgeData, 1, edgeNum, COPYOUT, fu); 
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_x0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, XU0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXU0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_x1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, XU1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXU1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_x2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, XU2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXU2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, mf);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, phiUDS);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, visac);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_D);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xpn0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXac0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xpn1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXac1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xpn2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, dPhidXac2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_n0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_n1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, face_n2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYOUT, rface_1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYOUT, rface_2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, rcpac);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, denac);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xac0);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xac1);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, Xac2);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, phiCDS);
	addSingleArray(frontEdgeData, 1, edgeNum, COPYIN, phiCDS2);
	constructEmptyArray(backEdgeData);
	constructEmptyArray(selfConnData);
	constructEmptyArray(vertexData);
	data.backEdgeData  = &backEdgeData;
	data.frontEdgeData = &frontEdgeData;
	data.selfConnData  = &selfConnData;
	data.vertexData    = &vertexData;
	cOpt[0].fun_slave = slave_swSeparate_vector;
	cOpt[0].fun_host  = swSeparate_vector;
	cOpt[0].data = &data;

	mlbIter->arrayIteration(&paraData, cOpt, 1);

	destroyArray(paraData);
	destroyArray(frontEdgeData);
	free(cOpt);
}

