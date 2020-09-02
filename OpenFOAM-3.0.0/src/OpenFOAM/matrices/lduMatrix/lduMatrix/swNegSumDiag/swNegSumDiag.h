#ifndef negSumDiag_H
#define negSumDiag_H
#include <stdlib.h>
#include <stdio.h>
#include "sw_struct.h"

#ifdef __cplusplus
extern "C" {
#endif
void negSumDiag_host(const SCALAR *lower,const SCALAR *upper,SCALAR *Diag,const LABEL *l,const LABEL *u,LABEL nCells,LABEL nFaces);
#ifdef __cplusplus
}
#endif

#endif
