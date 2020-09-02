#ifndef SEPARATE_VECTOR_HOST_H
#define SEPARATE_VECTOR_HOST_H
#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swSeparate_vector_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *mf, swFloat *phiUDS, swFloat *face_x0, 
			swFloat *dPhidXac0, swFloat *face_x1, swFloat *dPhidXac1,
			swFloat *face_x2, swFloat *dPhidXac2, swFloat *XU0, 
			swFloat *XU1, swFloat *XU2, swFloat *visac, swFloat *face_D,
			swFloat *Xpn0, swFloat *Xpn1, swFloat *Xpn2, swFloat *face_n0, 
			swFloat *face_n1, swFloat *face_n2, swFloat *dPhidXU0, 
			swFloat *dPhidXU1, swFloat *dPhidXU2, swFloat *fu, 
			swFloat *rface_1, swFloat *rface_2, swFloat *rcpac, 
			swFloat *denac, swFloat *Xac0, swFloat *Xac1, swFloat *Xac2, 
			swFloat *phiCDS, swFloat *phiCDS2, swFloat csBlend,
			swFloat vis_lam, swFloat vis_sigma, swFloat vis_PrScNr, 
			swFloat vis_lambda, swFloat vis_diffuse, swFloat vis_type,
			swFloat cvScheme, swFloat csFaceCorrect,
			swInt edgeNum, swInt vertexNum);

#endif
