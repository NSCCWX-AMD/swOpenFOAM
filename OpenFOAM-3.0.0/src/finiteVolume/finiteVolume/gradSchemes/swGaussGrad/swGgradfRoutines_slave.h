

#ifndef swGradfRoutines_slave_H
#define swGradfRoutines_slave_H

#include "rowSubsection.h"

#ifdef __cplusplus
extern "C"
{
#endif

struct swGgradf_scalar_paras
{
    const char* Sf_vptr; 
    const scalar* issf_ptr; 
    char* igGrad_vptr; 
    const label* lPtr; 
    const label* uPtr;
    const label fSize; 
    const label cSize; 
    const label vector_size; 
    const label vector_offset;
    const struct rowSubsection** secs;
    const label secNumInseg;
    const label colRoundNum;
};

struct swGgradf_vector_paras
{
    const char* Sf_vptr; 
    const char* issf_vptr; 
    char* igGrad_tptr; 
    const label* lPtr; 
    const label* uPtr;
    const label fSize; 
    const label cSize; 
    const label vector_size; 
    const label vector_offset;
    const label tensor_size;
    const label tensor_offset;
    const struct rowSubsection** secs;
    const label secNumInseg;
    const label colRoundNum;
};

typedef struct struct_swGgradfDivide_paras
{
    scalar* igGrad;
    const scalar* volume;
    const label NCell;
    const label size;
} swGgradfDivide_paras;

#ifdef __cplusplus
}
#endif


#endif
