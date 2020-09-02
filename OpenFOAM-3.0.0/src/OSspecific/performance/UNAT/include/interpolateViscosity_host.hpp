#ifndef INTERPOLATEVISCOSITY_HOST_HPP
#define INTERPOLATEVISCOSITY_HOST_HPP
#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swInterpolateViscosity_host(MultiLevelBlockIterator *mlbIter,
		   	swFloat *visSx, swFloat *visSy, swFloat *visSz, swFloat *visPNx,
			swFloat *visPNy, swFloat *visPNz, swFloat *rface0, 
			swFloat *rface1, swFloat *face_nx, swFloat *face_ny, 
			swFloat *face_nz, swFloat *face_dx, swFloat *face_dy, 
			swFloat *face_dz, swFloat *face_d, swFloat *viseff, 
			swFloat *facn, swInt edgeNum, swInt vertexNum);

#endif
