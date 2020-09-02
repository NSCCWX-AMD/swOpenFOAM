#include<stdio.h>
//#include "rowSubsection.h"
#ifndef SURFACEINTEGRATE_STRUCT_H
#define SURFACEINTEGRATE_STRUCT_H
 #define ArraySize 48000

#ifdef __cplusplus
extern "C"
{
#endif
  struct surfaceIntegrate_para {
            double  *ivf_Ptr;
            const double  *issf_Ptr;
            const int  *owner_Ptr;
            const int *neighbour_Ptr;
            const struct rowSubsection** secs;
            int secNumInSeg;
            int colRoundNum;
            int CEPs;
            int vector_size;
			int nFaces;
			int nCells;
          };


#ifdef __cplusplus
}
#endif

#endif
