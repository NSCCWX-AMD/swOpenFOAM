#ifndef USERFUNC_HOST_H
#define USERFUNC_HOST_H
#include "vectorOps_struct.h"

void vectorOps_host    		 (MVM_Arrays*, void (*userFuncPtr)(MVM_Arrays*));
void scalingFactor_host 	 (MVM_Arrays*);
void gSum_host				 (MVM_Arrays*, void (*userFuncPtr)(MVM_Arrays*));
void vectorCopy_host         (MVM_Arrays*);
void residualNormFactor_host (MVM_Arrays*);

#endif
