#include "slaveUtils.h"

void DMA_Get(void* dest, const void* source, const int size)
{
    DMA_Status status=0;
	if(size>0)
	{
	    dma_desc gv=0;
		A_DMA_GET_SET( gv, PE_MODE, size, &status);
	    A_DMA_GET_RUN( gv, source, dest);
		dma_wait(&status, 1);
	}
}

void DMA_IGet(void* dest, const void* source, const int size, 
        DMA_Status* status)
{
	if(size>0)
	{
	    dma_desc gv=0;
		A_DMA_GET_SET( gv, PE_MODE, size, status);
	    A_DMA_GET_RUN( gv, source, dest);
	}
}

void DMA_Put(void* dest, const void* source, const int size)
{
    DMA_Status status=0;
	if(size>0)
	{
		dma_desc pv=0;
	    A_DMA_PUT_SET( pv, PE_MODE, size, &status);
		A_DMA_PUT_RUN( pv, source, dest);
	    dma_wait(&status, 1);
	}
}

void DMA_IPut(void* dest, const void* source, const int size, 
        DMA_Status* status)
{
	if(size>0)
	{
	    dma_desc pv=0;
		A_DMA_PUT_SET( pv, PE_MODE, size, status);
	    A_DMA_PUT_RUN( pv, source, dest);
	}
}

void DMA_Wait( DMA_Status* status, const int num)
{
    dma_wait(status, num);
}

void athread_wait(int* reply, int count)
{
    while(*reply != count);
    return;
}
