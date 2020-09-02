#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include "surfaceIntegrate_struct.h"
#include"rowSubsection.h"

#include <dma.h>


#define ALIGNED(addr) ((((unsigned long)(addr)>>8)+1)<<8)

#define A_DMA_GET_SET(da,mode,len,re_addr)\
({  \
dma_set_op(&da,DMA_GET);    \
dma_set_mode(&da,mode); \
dma_set_size(&da,len);  \
dma_set_reply(&da,re_addr); \
})

#define A_DMA_GET_RUN(da,src,dest)  \
({  \
    dma(da,src,dest);   \
})

#define A_DMA_PUT_SET(da,mode,len,re_addr)\
({      \
dma_set_op(&da,DMA_PUT);        \
dma_set_mode(&da,mode); \
dma_set_size(&da,len);  \
dma_set_reply(&da,re_addr);     \
})

#define A_DMA_PUT_RUN(da,src,dest)      \
({      \
        dma(da,dest,src);       \
})

#define LDM_MALLOC(ptr,size) \
{ \
    if( get_allocatable_size() > size) \
        ptr = ldm_malloc(size); \
    else \
    { \
        printf("\n***Error: no enough LDM space to malloc\n"); \
        exit(-1); \
    } \
}

// can not be used in inbalance conditional branches
#define SERIAL_BEGIN \
{ \
    int serial_cpeid = athread_get_id(-1); \
    int serial_cpeIndex; \
    for(serial_cpeIndex = 0; \
        serial_cpeIndex < 64; \
        serial_cpeIndex++){\
        if(serial_cpeIndex == \
           serial_cpeid){

#define SERIAL_END \
        } \
        athread_syn(ARRAY_SCOPE, 0xFFFF); \
    } \
}

// can not be used in inbalance conditional branches
#define SPRINT(...) \
{ \
    int serial_cpeid = athread_get_id(-1); \
    int serial_cpeIndex; \
    for(serial_cpeIndex = 0; \
        serial_cpeIndex < 64; \
        serial_cpeIndex++){\
        if(serial_cpeIndex == \
           serial_cpeid){ \
            printf(__VA_ARGS__); \
        } \
        athread_syn(ARRAY_SCOPE, 0xFFFF); \
    } \
}

// can not be used in inbalance conditional branches
#define CPRINT(cpeID,...) \
{\
    if(athread_get_id(-1) == cpeID)\
        printf(__VA_ARGS__); \
    athread_syn(ARRAY_SCOPE, 0xFFFF);\
}


//__thread_local char data[ArraySize];
//__thread_local int  owner_slave[4000];
//__thread_local double  issf_slave[4000];
//__thread_local double  ivf_slave[4000];
//__thread_local int  my_id;

void surfaceIntegrate_slave(struct surfaceIntegrate_para *para)
{


    const double *issfPtr=para->issf_Ptr;
    double *ivfPtr=para->ivf_Ptr;
    const int  *ownerPtr=para->owner_Ptr;
    const int  *neighbourPtr=para->neighbour_Ptr;
    int nfaces=para->nFaces;
    int ncells=para->nCells;
    int vectorsize=para->vector_size;
    const struct rowSubsection** secs=para->secs;
    int  secNumInseg =para->secNumInSeg;
    int  colRoundNum= para->colRoundNum;

    int my_id=athread_get_id(-1);
    char data[ArraySize];

    int i;
    int j;
	int k;
	volatile int get_reply;
	volatile int put_reply;

    for(i=0;i<secNumInseg;i++)
    {
         const struct rowSubsection subsection = secs[my_id][i];
         int faceStart  = subsection.faceStart;
         int nFace= subsection.nFaces;
         char* dataEnd = data;

         volatile const int* owner = (int *)ALIGNED(dataEnd);
         dataEnd = owner + nFace;

        unsigned get_reply=0;
        dma_desc gv=0;
		A_DMA_GET_SET( gv, PE_MODE, nFace*sizeof(int), &get_reply);
        A_DMA_GET_RUN( gv, &ownerPtr[faceStart], owner);
        dma_wait(&get_reply, 1);

        int rowStart = owner[0];
        int rowNum   = owner[nFace-1] - owner[0] + 1;

	    const double* issf = (double *)ALIGNED(dataEnd);
        dataEnd = issf + nFace*vectorsize;

        double* ivf = (double *)ALIGNED(dataEnd);
        dataEnd = ivf + rowNum*vectorsize;

        if((int)( dataEnd - data) >= ArraySize)
        {
            printf("\n[%d]***Error: data size in ldm exceeds the bound in owner\n\n",
                    my_id);
            exit(1);
        }



		dma_desc gv1=0, gv2=0;
        volatile unsigned gReply1=0;
        A_DMA_GET_SET( gv1, PE_MODE, nFace*sizeof(double)*vectorsize, &gReply1);
		A_DMA_GET_RUN( gv1, &issfPtr[faceStart*vectorsize],issf);
        A_DMA_GET_SET( gv2, PE_MODE, rowNum*sizeof(double)*vectorsize, &gReply1);
        A_DMA_GET_RUN( gv2, &ivfPtr[rowStart*vectorsize], ivf);
        dma_wait(&gReply1, 2);


		 for(j=0;j<nFace;j++){
            for(k=0;k<vectorsize;k++){
          ivf[(owner[j]-owner[0])*vectorsize+k]+=issf[j*vectorsize+k];
            }
			}


	   dma_desc pv=0;
        volatile unsigned pReply=0;
        A_DMA_PUT_SET( pv, PE_MODE, rowNum*sizeof(double)*vectorsize, &pReply);
        A_DMA_PUT_RUN( pv, ivf, &ivfPtr[rowStart*vectorsize]);
        dma_wait(&pReply, 1);

    }

	 for(i=0; i<secNumInseg; i++)
    {
		const struct rowSubsection subsection = secs[my_id][i];
        const int faceStart  = subsection.faceStart;
        const int nFace      = subsection.nFaces;
        const int colRound   = subsection.colRound;
        const int nSecs      = subsection.nSecs;
        const int * colSAndC = subsection.colStartsAndCounts;

		int totalCols = 0;
        int colSecIndex = nSecs;
        while(colSecIndex != 0){
            totalCols += colSAndC[2*colSecIndex-1];
            colSecIndex--;
        }

		char* dataEnd = data;

        const double* issf = ALIGNED(dataEnd);
        dataEnd = issf + nFace*vectorsize;

        const int* neibs = ALIGNED(dataEnd);
        dataEnd = neibs + nFace;

		 dma_desc gv0=0, gv1=0, gv2=0;
        volatile unsigned  gReply1=0, gReply2=0;
        A_DMA_GET_SET( gv1, PE_MODE, nFace*sizeof(double)*vectorsize, &gReply1);
        A_DMA_GET_SET( gv2, PE_MODE, nFace*sizeof(int),  &gReply2);
        A_DMA_GET_RUN( gv1, &issfPtr[faceStart*vectorsize],            issf);
        A_DMA_GET_RUN( gv2, &neighbourPtr[faceStart],                neibs);

		double* tmp_ivf = ALIGNED(dataEnd);
        dataEnd = tmp_ivf + totalCols*vectorsize;

        double* ivf = ALIGNED(dataEnd);
        dataEnd = ivf + (unsigned) (0.4*nFace)*vectorsize;

		if((label)( dataEnd - data) >= ArraySize)
        {
            printf("\n[%d]***Error: data size in ldm exceeds the bound in neighbour\n\n",
                    my_id);
            exit(1);
        }

		  memset(tmp_ivf, 0, totalCols*vectorsize*sizeof(double));

        dma_wait(&gReply1, 1);
        dma_wait(&gReply2, 1);

		 int faceIndex;
        for(faceIndex=0; faceIndex<nFace; faceIndex++){
          //  const scalar * vSf = &Sf[ (faceIndex)*vectorsize + vector_offset ];

            const int col_global = neibs[faceIndex];
            int col_local = 0;
            int secI = 0;
            while( col_global < colSAndC[2*secI] ||
                   col_global >= colSAndC[2*secI+1]+colSAndC[2*secI]
                 )
            {
                if( secI < nSecs ){
                    col_local += colSAndC[2*secI+1];
                    secI++;
                }
                else{
                    printf("[%d]***Error: col_global %d, original faceIndex %d, original col_global %d\n",
                            my_id, neibs[faceIndex], faceStart+faceIndex, neighbourPtr[faceStart+faceIndex]);
                    exit(-1);
                }
            }

            col_local += col_global - colSAndC[2*secI];
           for(k=0;k<vectorsize;k++)
            tmp_ivf[col_local*vectorsize + k] -= issf[faceIndex*vectorsize + k];

            //scalar * vigGrad = &igGrad_vptr[ neibs[faceIndex]*vector_size + vector_offset ];
            //vigGrad[0] += vSf[0]*issf[ faceIndex ];
            //vigGrad[1] += vSf[1]*issf[ faceIndex ];
            //vigGrad[2] += vSf[2]*issf[ faceIndex ];
        }

		int round;
        for(round = 0; round < colRoundNum; round++){
            athread_syn( ARRAY_SCOPE, 0xFFFF);
            if( round == colRound )
            {
                int secI;
                int tmpStart = 0;
                for( secI = 0; secI < nSecs; secI++){
                    int colStart = colSAndC[2*secI];
                    int colNum   = colSAndC[2*secI+1];
					/*
					if(my_id==0)
					{
						printf("%d %d\n",colStart,colNum);

					}*/

                    dma_desc gv3=0;
                    volatile unsigned gReply3=0;
                    A_DMA_GET_SET( gv3, PE_MODE, colNum*vectorsize*sizeof(double), &gReply3);
                    A_DMA_GET_RUN( gv3, &ivfPtr[(colStart)*vectorsize], ivf);
                    dma_wait(&gReply3, 1);

					/*int t;
					if(my_id==0)
						for(t=0;t<colNum*vectorsize;t++){
							printf("%f",ivf[t]);
							if(t%10==0)
								printf("\n");
						}*/

						//return;

                    int cellIndex;
                    for(cellIndex=0; cellIndex<colNum; cellIndex++)
                    {
                        double* vc = &ivf[cellIndex*vectorsize];
						for(k=0;k<vectorsize;k++)
                        vc[k] += tmp_ivf[ (cellIndex+tmpStart)*vectorsize + k];

                    }


                   dma_desc pv=0;
                    volatile unsigned pReply=0;
                    A_DMA_PUT_SET( pv, PE_MODE, colNum*vectorsize*sizeof(double), &pReply);
                    A_DMA_PUT_RUN( pv, ivf, &ivfPtr[colStart*vectorsize]);
                    dma_wait(&pReply, 1);



                    tmpStart += colNum;
                }
            }
        }



	}




}





