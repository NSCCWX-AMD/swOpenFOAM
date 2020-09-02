#ifndef SEPARATE_INTERPOLATION_HOST_HPP
#define SEPARATE_INTERPOLATION_HOST_HPP
#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT; 

void swSeparate_interpolation_host(MultiLevelBlockIterator *mlbIter,
			swFloat* face_lam, swFloat *vis,
			swFloat *visac, swFloat *gradPhi0, swFloat *dPhidXac0, 
			swFloat *gradPhi1, swFloat *dPhidXac1, swFloat *gradPhi2,
			swFloat *dPhidXac2, swFloat *mf, swFloat *cell_x0, 
			swFloat *XU0, swFloat *cell_x1, swFloat *XU1, swFloat* cell_x2,
		   	swFloat *XU2, swFloat *dPhidXU0, swFloat *dPhidXU1, 
			swFloat *dPhidXU2, swFloat *phi, swFloat *phiUDS, swFloat *Xpn0,
		   	swFloat *Xpn1, swFloat *Xpn2, swFloat *phiCDS, swFloat *Xac0, 
			swFloat *Xac1, swFloat *Xac2, swFloat *rcp, swFloat *rcpac, 
			swFloat *den, swFloat *denac, swFloat *phiCDS2, 
			swInt edgeNum, swInt vertexNum);

#endif
