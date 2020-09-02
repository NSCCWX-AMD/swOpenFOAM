#ifndef MULTILEVELBLOCKITERATOR_H
#define MULTILEVELBLOCKITERATOR_H

#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#include "RlmpiInitializer.hxx"
#endif

typedef struct 
{
	// topology relevent
	swInt* owner;
	swInt* neighbor;
	swInt* cellStarts;
	swInt* blockStarts;
	swInt* recvStarts;
	swInt  vertexNumber;
	swInt  edgeNumber;
	swInt  mshBlockNum;
	swInt  maxEdges;
	swInt  maxCells;

	// Run-time data
	swInt spIndex;

	// Register communication topology
	RlmpiInfo *schedule_data;

	// computing data
	Arrays* paraData;

	coupledOperator *cOpt;
	int optNum;

} MLB_edge2VertexPara;

typedef struct 
{
	// topology relevent
	swInt* firstEdgeVertice;
	swInt* vertexNeighbor;
	swInt* cellStarts;
	swInt* blockStarts;
	swInt  vertexNumber;
	swInt  edgeNumber;
	swInt  mshBlockNum;
	swInt  maxEdges;
	swInt  maxCells;

	// Run-time data
	swInt spIndex;

	// Register communication topology
	RlmpiInfo *schedule_data;

	// computing data
	Arrays* paraData;

	coupledOperator *cOpt;
	int optNum;

} MLB_vertex2EdgePara;

typedef struct 
{
	// topology relevent
	swInt* cellStarts;
	swInt* blockStarts;
	swInt  vertexNumber;
	swInt  edgeNumber;
	swInt  mshBlockNum;

	// Run-time data
	swInt spIndex;

	// computing data
	Arrays* paraData;

	coupledArrayOperator *cOpt;
	int optNum;

} MLB_arrayPara;

#ifdef __cplusplus
}
#endif



#endif // MULTILEVELBLOCKITERATOR_H
