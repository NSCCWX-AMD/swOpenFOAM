#ifndef ROWSUBSECTIONITERATOR_H
#define ROWSUBSECTIONITERATOR_H

#include "iterator.h"
#include "rowSubsection.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct 
{
	// topology relevent
	swInt* firstEdgeVertice;
	swInt* vertexNeighbor;
	const struct rowSubsection** secs;
	swInt  secNumInSeg;
	swInt colRoundNum;
	swInt  vertexNumber;
	swInt  edgeNumber;
	swInt* maxFacesInSeg;
	swInt* maxCellsInSeg;
	swInt* maxColsInSeg;
	swInt  maxColNum;

	// computing data
	Arrays* paraData;

	coupledOperator *cOpt;
	int optNum;

} RSS_vertex2EdgePara;


#ifdef __cplusplus
}
#endif



#endif // ROWSUBSECTIONITERATOR_H
