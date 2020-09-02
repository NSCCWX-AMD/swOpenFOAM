#ifndef INTEGRATE_H
#define INTEGRATE_H
#include "swMacro.h"
#include "iterator.h"
#include <sys/time.h>

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(swSurfaceIntegrate);
define_e2v_slaveFunPtr(swSurfaceIntegrate);

#define swSurfaceIntegrate_host(owner, neighbor, issf, ivf, \
			edgeNum, vertexNum, dims) \
{ \
	double time1,time2; \
	getTime(time1); \
	static UNAT::Topology* topo = UNAT::Topology::constructFromEdge( \
				owner, neighbor, edgeNum); \
	static std::vector<swInt> cellWeights(topo->getVertexNumber(), 2*dims);\
	static std::vector<swInt> edgeWeights(topo->getEdgeNumber(), 2*dims); \
	static UNAT::DirectSegmentIterator iterator(*topo, \
				&cellWeights[0], &edgeWeights[0]); \
	\
	static Arrays backEdgeData, frontEdgeData, selfConnData, vertexData; \
	constructSingleArray(backEdgeData, dims,edgeNum,COPYIN,issf);\
	constructSingleArray(frontEdgeData,dims,edgeNum,COPYIN,issf);\
	constructEmptyArray(selfConnData);\
	constructSingleArray(vertexData,dims,vertexNum,COPYOUT,ivf);\
	iterator.edge2VertexIteration(&backEdgeData,&frontEdgeData,\
				&selfConnData,&vertexData, \
				swSurfaceIntegrate,slave_swSurfaceIntegrate);\
	getTime(time2); \
} \

//	printf("slave core Compute: %f us\n", (time2-time1)*1000000); 
//void slave_integrate(Arrays* backEdgeData, Arrays* frontEdgeData,
//			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
//			swInt* endVertices);
//void integrate(Arrays* backEdgeData, Arrays* frontEdgeData,
//			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
//			swInt* endVertices);
#ifdef __cplusplus
}
#endif

#endif
