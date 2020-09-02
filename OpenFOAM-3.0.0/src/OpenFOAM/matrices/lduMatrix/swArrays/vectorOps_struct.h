#ifndef vectorOps_struct_H
#define vectorOps_struct_H

#define SCALAR double
#define LABEL int

#define ALIGNED(addr) ((( ( (unsigned long)(addr) - 1)>>5)+1)<<5)
#define ArraySize 57344

#define slaveDoubleNumber 1536
#define accUsingSize 20000

typedef struct
{
    SCALAR *A1Ptr;      // default set as result, keep it NULL if there is no result return
    SCALAR *A2Ptr;
    SCALAR *A3Ptr;
    SCALAR *A4Ptr;
    SCALAR *k1Ptr;
    SCALAR *k2Ptr;
    SCALAR  k1;         // parameter 1
    SCALAR  k2;         // parameter 2
    LABEL   n1;         // size of array
    LABEL   nCells;     // size of local loop
    LABEL   returnA2;   // if return A2
}MVM_Arrays;

typedef struct
{
    SCALAR *A1Ptr;      // default set as result
    SCALAR *A2Ptr;
    SCALAR *A3Ptr;
    SCALAR *A4Ptr;
    SCALAR *k1Ptr;      // for scaling factor
    SCALAR *k2Ptr;      // for scaling factor
    SCALAR  k1;         // parameter 1
    SCALAR  k2;         // parameter 2
    LABEL   n1;         // size of array
    LABEL   returnA2;   // if return A2
    void  (*userFuncPtr)(MVM_Arrays*);
}MVM_ParametersPack;

inline void init_MVM_Arrays(MVM_Arrays *MVM_ArraysPtr, LABEL size)
{
    MVM_ArraysPtr->A1Ptr    = 0;
    MVM_ArraysPtr->A2Ptr    = 0;
    MVM_ArraysPtr->A3Ptr    = 0;
    MVM_ArraysPtr->A4Ptr    = 0;
    MVM_ArraysPtr->k1Ptr    = 0;
    MVM_ArraysPtr->k2Ptr    = 0;
    MVM_ArraysPtr->k1       = 0.0;
    MVM_ArraysPtr->k2       = 0.0;
    MVM_ArraysPtr->n1       = size;
    MVM_ArraysPtr->nCells   = 0;
    MVM_ArraysPtr->returnA2 = 0;
}

#endif
