#include <slave.h>
#include <dma.h>

#ifndef SLAVEUTILS_H
#define SLAVEUTILS_H

#define ALIGNED(addr) ((((unsigned long)(addr-1)>>5)+1)<<5)

#define A_DMA_GET_SET(da,mode,len,re_addr) \
{ \
    dma_set_op(&da,DMA_GET); \
    dma_set_mode(&da,mode); \
    dma_set_size(&da,len); \
    dma_set_reply(&da,re_addr); \
}

#define A_DMA_GET_RUN(da,src,dest) \
{ \
    dma(da,src,dest); \
}

#define A_DMA_PUT_SET(da,mode,len,re_addr) \
{ \
    dma_set_op(&da,DMA_PUT); \
    dma_set_mode(&da,mode); \
    dma_set_size(&da,len); \
    dma_set_reply(&da,re_addr); \
}

#define A_DMA_PUT_RUN(da,src,dest) \
{ \
    dma(da,dest,src); \
}


typedef unsigned long DMA_Status;

void DMA_Get(void* dest, const void* source, const int size);

void DMA_IGet(void* dest, const void* source, const int size, 
        DMA_Status* status);

void DMA_Put(void* dest, const void* source, const int size);

void DMA_IPut(void* dest, const void* source, const int size, 
        DMA_Status* status);

void DMA_Wait( DMA_Status* status, const int num);

void athread_wait(int* reply, int count);

#endif
