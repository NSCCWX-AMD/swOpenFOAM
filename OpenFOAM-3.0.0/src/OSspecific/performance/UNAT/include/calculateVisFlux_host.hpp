#ifndef CALCULATEVISFLUX_HOST_HPP
#define CALCULATEVISFLUX_HOST_HPP

#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swCalculateVisFlux_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *facn, swFloat *visSy, swFloat *dudy, swFloat *dvdx, 
			swFloat *dvdy, swFloat *dvdz, swFloat *dwdy, swFloat *Su, 
			swFloat *Sv, swFloat *Sw, swFloat *dudx, swFloat *dudz, 
			swFloat *dwdx, swFloat *dwdz, swFloat *visSx, swFloat *visSz,
			swFloat *visPNx, swFloat *visPNy, swFloat *visPNz, 
			swInt edgeNum, swInt vertexNum);

#endif
