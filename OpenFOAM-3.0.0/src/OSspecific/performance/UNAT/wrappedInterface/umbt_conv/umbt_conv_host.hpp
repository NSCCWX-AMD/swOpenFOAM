#ifndef UMBT_CONV_HOST_HPP
#define UMBT_CONV_HOST_HPP

#include "swMacro.h"
#include "iterator.hpp"

using namespace UNAT;

void umbt_conv_host(Iterator *iter, swFloat *n_vg, swFloat *u, swFloat *cen,
			swFloat *gra, swFloat *rhs, swInt edgeNum, swInt cellNum);

#endif
