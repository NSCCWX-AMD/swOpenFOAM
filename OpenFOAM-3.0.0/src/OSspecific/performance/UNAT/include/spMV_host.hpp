#ifndef SPMV_HOST_HPP
#define SPMV_HOST_HPP

#include "swMacro.h"
#include "iterator.hpp"

using namespace UNAT;

void swSpMV_v2e_host(Iterator *iter, swFloat *edge, 
			swFloat *bPtr, swFloat *xPtr, swFloat *diagPtr, 
			swInt dims, swInt edgeNum, swInt cellNum);
void swSpMV_e2v_host(Iterator *iter, swFloat *lowerPtr, swFloat *upperPtr, 
			swFloat *bPtr, swFloat *xPtr, swFloat *diagPtr,
			swInt dims, swInt edgeNum, swInt cellNum);

#endif
