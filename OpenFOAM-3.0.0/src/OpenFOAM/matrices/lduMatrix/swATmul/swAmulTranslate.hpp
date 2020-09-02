#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "amulMacros.h"
extern "C" {
}

refilltion translate_matrix_matrix(amul_translate_array_ptr &data_array,
								   amul_para_ptr matrix,
								   LABEL size,
								   const LABEL offset_size);


void compare_debug(amul_translate_array_ptr data_array,
				   const LABEL *uPtr,
				   const LABEL *lPtr,
				   const double *upperPtr,
				   const double *lowerPtr,
				   const LABEL nFaces,
				   LABEL offset_size);

void translate_refill_data_host(refilltion_ptr parameter);
amul_translate_array_ptr translate_matrix_matrix_T(const amul_translate_array_ptr data_array);

#define translate_refill_data(para1) translate_refill_data_host(para1)

