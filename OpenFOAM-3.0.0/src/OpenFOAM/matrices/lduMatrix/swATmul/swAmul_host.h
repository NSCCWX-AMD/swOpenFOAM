#ifndef SW_AMUL_HOST_H
#define SW_AMUL_HOST_H
#include "amulMacros.h"

void amul_host(amul_para_ptr in,amul_translate_array_ptr in2);
#define tmul_host(para1,para2) amul_host(para1,para2)


#endif
