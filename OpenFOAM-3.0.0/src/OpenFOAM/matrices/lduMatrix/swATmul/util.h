#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include "amulMacros.h"
LABEL binarySearch(const LABEL *arr, LABEL start, LABEL end, LABEL target,LABEL offset_size) ;//supremum

void get_block_row_end_index(const LABEL * arr , LABEL * res , LABEL size , LABEL block_row_size,LABEL offset_size) ;

#define SIMPLE_MALLOC(type,size) ((type*)malloc(sizeof(type)*(size)))
#endif
