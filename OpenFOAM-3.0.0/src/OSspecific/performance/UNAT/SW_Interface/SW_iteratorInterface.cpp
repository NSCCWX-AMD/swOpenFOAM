#include "SW_iteratorInterface.h"

#include "topology.hpp"
#include "iterator.hpp"
#include "multiLevelBlockIterator.hpp"
#include "rowSubsectionIterator.hpp"
// using namespace UNAT;

// UNAT::Iterator *iter;

void constructMLBIterator(
			swInt* rowAddr, swInt* colAddr, 
			swInt* cellWeights, swInt* faceWeights,
			swInt faceNum)
{
	// UNAT::Topology *topo = UNAT::Topology::constructFromEdge(rowAddr, colAddr, faceNum);
	// iter = new UNAT::MultiLevelBlockIterator(*topo, cellWeights,
		// faceWeights, true);
}


// swInt* getEdgeMap(UNAT::MultiLevelBlockIterator *mlbIter)
// {
	// return (mlbIter->getPostEdgeOrder());
// }
// 
// swInt* getCellMap(UNAT::MultiLevelBlockIterator *mlbIter)
// {
	// return (mlbIter->getPostVertexOrder());
// }

void constructRSSIterator(
			swInt* rowAddr, swInt* colAddr, 
			swInt* cellWeights, swInt* faceWeights,
			swInt  faceNum, Arrays* paraData, coupledOperator* opt,
			swInt  optNum)
{

	// printf("constructRSSIterator\n");
	// Topology *topo;
	// int test = topo->getVertexNumber();
	static UNAT::Topology *topo = UNAT::Topology::constructFromEdge(rowAddr, colAddr, faceNum);
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(opt, *topo, cellWeights,
		faceWeights, E2V, true);
	iter->edge2VertexIteration(paraData, opt, optNum);
}