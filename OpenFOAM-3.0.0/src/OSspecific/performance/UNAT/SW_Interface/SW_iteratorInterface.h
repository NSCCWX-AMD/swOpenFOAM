#ifndef SW_ITERATORINTERFACE_H
#define SW_ITERATORINTERFACE_H

#include "swMacro.h"
// #include "register.H"
#include "iterator.h"



// swInt* getEdgeMap(UNAT::MultiLevelBlockIterator *mlbIter);
// swInt* getCellMap(UNAT::MultiLevelBlockIterator *mlbIter);

#ifdef __cplusplus
extern "C"
{
#endif
void constructRSSIterator(
			swInt* rowAddr, swInt* colAddr,
			swInt* cellWeights, swInt* faceWeights,
			swInt  faceNum, Arrays* paraData, coupledOperator* opt,
			swInt  optNum);

void constructMLBIterator(
			swInt* rowAddr, swInt* colAddr, 
			swInt* cellWeights, swInt* faceWeights,
			swInt  faceNum);
#ifdef __cplusplus
}
#endif



#endif
