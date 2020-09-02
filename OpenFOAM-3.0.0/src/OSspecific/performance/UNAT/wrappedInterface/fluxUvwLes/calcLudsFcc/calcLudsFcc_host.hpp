#ifndef CALCLUDSFCC_HOST_HPP
#define CACLLUDSFCC_HOST_HPP
#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swCalcLudsFcc_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *massFlux, swFloat *facex_x, swFloat *facex_y, 
			swFloat *facex_z, swFloat *cellx_x, swFloat *cellx_y, 
			swFloat *cellx_z, swFloat *fccx, swFloat *fccy, swFloat *fccz,
			swFloat *rface0, swFloat *rface1, swFloat *facn, swFloat CS, 
			swFloat gamblend, swInt edgeNum, swInt vertexNum);

#endif
