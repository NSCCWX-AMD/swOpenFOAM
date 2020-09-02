#ifndef CALCULATEFCCFLUX_HOST_HPP
#define CALCULATEFCCFLUX_HOST_HPP

#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swCalculateFccFlux_host(MultiLevelBlockIterator *mlbIter, 
		   	swFloat *massFlux, swFloat *fccx, swFloat *dudx, swFloat *dvdx,
		   	swFloat *dwdx, swFloat *Su, swFloat *Sv, swFloat *Sw, 
			swFloat *fccy, swFloat *fccz, swFloat *dudy, swFloat *dvdy, 
			swFloat *dwdy, swFloat *dudz, swFloat *dvdz, swFloat *dwdz, 
			swInt edgeNum, swInt vertexNum);

#endif
