#ifndef SEPARATE_INTEGRATE_HOST_HPP
#define SEPARATE_INTEGRATE_HOST_HPP
#include "swMacro.h"
#include "iterator.h"
#include "multiLevelBlockIterator.hpp"

using namespace UNAT;

void swSeparateIntegrate_host(MultiLevelBlockIterator *mlbIter, 
			swFloat *fu, swFloat *su, swInt edgeNum, swInt vertexNum);

#endif
