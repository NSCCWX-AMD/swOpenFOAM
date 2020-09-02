

#ifndef swGgradfRoutines_H
#define swGgradfRoutines_H

#include "rowSubsection.h"

#ifdef __cplusplus
extern "C"
{
#endif

void swGgradfInnerRoutine_scalar( 
        const char* Sf_vptr, 
        const scalar* issf_ptr, 
        char* igGrad_vptr, 
        const label* lPtr, 
        const label* uPtr,
        const label fSize, 
        const label cSize, 
        const label vector_size, 
        const label vector_offset,
        const struct rowSubsection** secs,
        const label secNumInseg,
        const label colRoundNum
        );


void swGgradfInnerRoutine_vector( 
        const char* Sf_vptr, 
        const char* issf_vptr, 
        char* igGrad_tptr, 
        const label* lPtr, 
        const label* uPtr,
        const label fSize, 
        const label cSize, 
        const label vector_size, 
        const label vector_offset,
        const label tensor_size, 
        const label tensor_offset,
        const struct rowSubsection** secs,
        const label secNumInseg,
        const label colRoundNum
        );

void swGgradfDivide_host(
        scalar* igGrad,
        const scalar* volume,
        const label NCell,
        const label size
        );
#ifdef __cplusplus
}
#endif

#endif
