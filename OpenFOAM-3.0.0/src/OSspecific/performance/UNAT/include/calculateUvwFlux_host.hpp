#ifndef CALCULATEUVWFLUX_HOST_HPP
#define CALCULATEUVWFLUX_HOST_HPP
#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swCalculateUvwFlux_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *massFlux, swFloat *u, swFloat *v, swFloat *w, 
			swFloat *Su, swFloat *Sv, swFloat *Sw, swFloat gamblend, 
			swInt edgeNum, swInt vertexNum);

#endif
