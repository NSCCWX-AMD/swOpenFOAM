#ifndef NEGSUMDIAG_H
#define NEGSUMDIAG_H
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(swNegSumDiag);
define_e2v_slaveFunPtr(swNegSumDiag);

#define swNegSumDiag_host(owner, neighbor, lowerPtr, upperPtr, diagPtr, \
			edgeNum, vertexNum)  \
{ \
	getTime(time1);\
	static UNAT::Topology* topo = UNAT::Topology::constructFromEdge(owner,neighbor,edgeNum);\
	static std::vector<swInt> cellWeights(topo->getVertexNumber(), 2);\
	static std::vector<swInt> edgeWeights(topo->getEdgeNumber(), 2);\
	static UNAT::DirectSegmentIterator iterator(*topo, &cellWeights[0], &edgeWeights[0]);\
	\
	static Arrays backEdgeData, frontEdgeData, selfConnData, vertexData;\
	constructSingleArray(backEdgeData,1,edgeNum,COPYIN,upperPtr);\
	constructSingleArray(frontEdgeData,1,edgeNum,COPYIN,lowerPtr);\
	constructEmptyArray(selfConnData);\
	constructSingleArray(vertexData,1,vertexNum,COPYOUT,diagPtr);\
	iterator.edge2VertexIteration(&backEdgeData,&frontEdgeData,\
				&selfConnData,&vertexData,swNegSumDiag,slave_swNegSumDiag);\
	getTime(time2);\
} \

//printf("slave core Compute: %f us\n", (time2-time1)*1000000);
#ifdef __cplusplus
}
#endif

#endif
