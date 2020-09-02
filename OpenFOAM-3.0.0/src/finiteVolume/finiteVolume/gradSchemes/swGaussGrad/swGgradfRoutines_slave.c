#include "swGgradfRoutines_slave.h"
#include "stdlib.h"
#include <slave.h>
#include <dma.h>
#include "slaveUtils.h"

#define ALIGNED(addr) ((((unsigned long)(addr-1)>>5)+1)<<5)

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
}

// declare the global exception signal array
extern int CG_signals[];

//__thread_local char* data;

int swGgradfInner_scalar_slave(struct swGgradf_scalar_paras * _paras)
{
	char data[48000];

    struct swGgradf_scalar_paras paras = *_paras;
    const char* Sf_vptr               = paras.Sf_vptr;
    const scalar* issf_ptr            = paras.issf_ptr;
    char* igGrad_vptr                 = paras.igGrad_vptr;
    const label* lPtr                 = paras.lPtr;
    const label* uPtr                 = paras.uPtr;
    const label fSize                 = paras.fSize;
    const label cSize                 = paras.cSize;
    const label vector_size           = paras.vector_size;
    const label vector_offset         = paras.vector_offset;
    const struct rowSubsection** secs = paras.secs;
    const label secNumInseg           = paras.secNumInseg;
    const label colRoundNum           = paras.colRoundNum;

    const label cpeId = athread_get_id(-1);
    //__thread_local char data[48000];

    label secIndex;
    for(secIndex=0; secIndex< secNumInseg; secIndex++)
    {

        const struct rowSubsection subsection = secs[cpeId][secIndex];
        const label faceStart  = subsection.faceStart;
        const label nFace      = subsection.nFaces;

        /* get the row index and calculate row start and number */
        char* dataEnd = (char*)data;

        const label* owner = (const label*)ALIGNED(dataEnd);
        dataEnd = owner + nFace;

        volatile unsigned long get_reply=0;
        dma_desc gv=0;
        A_DMA_GET_SET( gv, PE_MODE, nFace*sizeof(label), &get_reply);
        A_DMA_GET_RUN( gv, &lPtr[faceStart], owner);
        DMA_Wait(&get_reply, 1);

        const label rowStart = owner[0];
        const label rowNum   = owner[nFace-1] - owner[0] + 1;

        // ldm_malloc memory
        const char* Sf = (const char*)ALIGNED(dataEnd);
        dataEnd = Sf + nFace*vector_size;

        const scalar* issf = (const scalar*)ALIGNED(dataEnd);
        dataEnd = issf + nFace;

        char* igGrad = (char*)ALIGNED(dataEnd);
        dataEnd = igGrad + rowNum*vector_size;

        if((label)( dataEnd - data) >= 48000)
        {
            CG_signals[cpeId] = 1;
            // core group must syncronize for colRoundNum
            // times in case normal cores block after bad
            // cores return
            label round;
            for(round = 0; round < colRoundNum*secNumInseg; round++)
                athread_syn( ARRAY_SCOPE, 0xFFFF);
            return 0;
        }

        // get the face fields
        dma_desc gv0=0, gv1=0, gv2=0;
        volatile unsigned long gReply0=0, gReply1=0, gReply2=0;

        A_DMA_GET_SET( gv0, PE_MODE, nFace*vector_size,    &gReply0);
        A_DMA_GET_SET( gv1, PE_MODE, nFace*sizeof(scalar), &gReply1);
        A_DMA_GET_SET( gv2, PE_MODE, rowNum*vector_size, &gReply2);
        A_DMA_GET_RUN( gv0, &Sf_vptr[faceStart*vector_size],    Sf);
        A_DMA_GET_RUN( gv1, &issf_ptr[faceStart],               issf);
        A_DMA_GET_RUN( gv2, &igGrad_vptr[rowStart*vector_size], igGrad);

        DMA_Wait(&gReply0, 1);
        DMA_Wait(&gReply1, 1);
        DMA_Wait(&gReply2, 1);

        // calculate the tmp grad
        label faceIndex;
        for(faceIndex=0; faceIndex<nFace; faceIndex++){
            scalar* vSf     = &Sf[ faceIndex*vector_size + vector_offset ];
            scalar* vigGrad = &igGrad[ (owner[faceIndex]-rowStart)*vector_size + vector_offset ];

            vigGrad[0] += vSf[0]*issf[faceIndex];
            vigGrad[1] += vSf[1]*issf[faceIndex];
            vigGrad[2] += vSf[2]*issf[faceIndex];
        }

        // add the tmp igGrad back
        dma_desc pv=0;
        volatile unsigned long pReply=0;
        A_DMA_PUT_SET( pv, PE_MODE, rowNum*vector_size, &pReply);
        A_DMA_PUT_RUN( pv, igGrad, &igGrad_vptr[rowStart*vector_size]);
        DMA_Wait(&pReply, 1);

    }
    /* calculate the column projection */
    for(secIndex=0; secIndex<secNumInseg; secIndex++)
    {

        const struct rowSubsection subsection = secs[cpeId][secIndex];
        const label faceStart  = subsection.faceStart;
        const label nFace      = subsection.nFaces;
        const label colRound   = subsection.colRound;
        const label nSecs      = subsection.nSecs;
//        const label * colSAndC = subsection.colStartsAndCounts;
		const label colSAndC[nSecs*2];
		dma_desc gv=0;
		volatile unsigned long gReply=0;
		A_DMA_GET_SET(gv,PE_MODE,nSecs*sizeof(label)*2,&gReply);
		A_DMA_GET_RUN(gv,subsection.colStartsAndCounts,&colSAndC[0]);
		DMA_Wait(&gReply,1);

        label totalCols = 0;
        label colSecIndex = nSecs;
        while(colSecIndex != 0){
            totalCols += colSAndC[2*colSecIndex-1];
            colSecIndex--;
        }

        char* dataEnd = (char*)data;

        const char* Sf = (const char*)ALIGNED(dataEnd);
        dataEnd = Sf + nFace*vector_size;

        const scalar* issf = (const scalar*)ALIGNED(dataEnd);
        dataEnd = issf + nFace;

        const label* neibs = (const label*)ALIGNED(dataEnd);
        dataEnd = neibs + nFace;

        dma_desc gv0=0, gv1=0, gv2=0;
        volatile unsigned long gReply0=0, gReply1=0, gReply2=0;

        A_DMA_GET_SET( gv0, PE_MODE, nFace*vector_size,    &gReply0);
        A_DMA_GET_SET( gv1, PE_MODE, nFace*sizeof(scalar), &gReply1);
        A_DMA_GET_SET( gv2, PE_MODE, nFace*sizeof(label),  &gReply2);
        A_DMA_GET_RUN( gv0, &Sf_vptr[faceStart*vector_size], Sf);
        A_DMA_GET_RUN( gv1, &issf_ptr[faceStart],            issf);
        A_DMA_GET_RUN( gv2, &uPtr[faceStart],                neibs);

        scalar* tmp_igGrad = (scalar*)ALIGNED(dataEnd);
        dataEnd = tmp_igGrad + totalCols*3;

        char* igGrad = (char*)ALIGNED(dataEnd);
        dataEnd = igGrad + 1024; //(unsigned) (0.4*nFace)*vector_size;
        if((label)( dataEnd - data) >= 48000)
        {
            CG_signals[cpeId] = 2;
            // core group must syncronize for colRoundNum
            // times in case normal cores block after bad
            // cores return
            label round;
            for(
                    round = 0;
                    round < colRoundNum*(secNumInseg-secIndex);
                    round++
                )
                athread_syn( ARRAY_SCOPE, 0xFFFF);
            return 0;
        }
        memset(tmp_igGrad, 0, totalCols*3*sizeof(scalar));

        DMA_Wait(&gReply0, 1);
        DMA_Wait(&gReply1, 1);
        DMA_Wait(&gReply2, 1);

        // calculate the tmp grad
        label faceIndex;
        for(faceIndex=0; faceIndex<nFace; faceIndex++){
            const scalar * vSf = &Sf[ (faceIndex)*vector_size + vector_offset ];

            const label col_global = neibs[faceIndex];
            label col_local = 0;
            label secI = 0;
            while( col_global < colSAndC[2*secI] ||
                   col_global >= colSAndC[2*secI+1]+colSAndC[2*secI]
                 )
            {
                if( secI < nSecs ){
                    col_local += colSAndC[2*secI+1];
                    secI++;
                }
                else{
                        CG_signals[cpeId] = 3;
                        // core group must syncronize for colRoundNum
                        // times in case normal cores block after bad
                        // cores return
                        label round;
                        for(
                                round = 0;
                                round < colRoundNum*(secNumInseg-secIndex);
                                round++
                            )
                            athread_syn( ARRAY_SCOPE, 0xFFFF);
                        return 0;
                }
            }

            col_local += col_global - colSAndC[2*secI];

            tmp_igGrad[col_local*3 + 0] -= vSf[0]*issf[faceIndex];
            tmp_igGrad[col_local*3 + 1] -= vSf[1]*issf[faceIndex];
            tmp_igGrad[col_local*3 + 2] -= vSf[2]*issf[faceIndex];

        }

        // communicate and write back according colRound
        label round;
        for(round = 0; round < colRoundNum; round++){
            athread_syn( ARRAY_SCOPE, 0xFFFF);
            if( round == colRound )
            {
                label secI;
                label tmpStart = 0;
                for( secI = 0; secI < nSecs; secI++){
                    label colStart = colSAndC[2*secI];
                    label colNum   = colSAndC[2*secI+1];

                    dma_desc gv3=0;
                    volatile unsigned long gReply3=0;
                    A_DMA_GET_SET( gv3, PE_MODE, colNum*vector_size, &gReply3);
                    A_DMA_GET_RUN( gv3, &igGrad_vptr[colStart*vector_size], igGrad);
                    DMA_Wait(&gReply3, 1);

                    label cellIndex;
                    for(cellIndex=0; cellIndex<colNum; cellIndex++)
                    {
                        scalar* vc = &igGrad[cellIndex*vector_size + vector_offset];
                        vc[0] += tmp_igGrad[ (cellIndex+tmpStart)*3 + 0];
                        vc[1] += tmp_igGrad[ (cellIndex+tmpStart)*3 + 1];
                        vc[2] += tmp_igGrad[ (cellIndex+tmpStart)*3 + 2];
                    }

                    dma_desc pv=0;
                    volatile unsigned long pReply=0;
                    A_DMA_PUT_SET( pv, PE_MODE, colNum*vector_size, &pReply);
                    A_DMA_PUT_RUN( pv, igGrad, &igGrad_vptr[colStart*vector_size]);
                    DMA_Wait(&pReply, 1);

                    tmpStart += colNum;
                }
            }
        }
    }
    return 0;
}




int swGgradfInner_vector_slave(struct swGgradf_vector_paras* _paras)
{
	char data[48000];

    struct swGgradf_vector_paras paras = *_paras;
    const char* Sf_vptr                 = paras.Sf_vptr;
    const char* issf_vptr               = paras.issf_vptr;
    char* igGrad_tptr                   = paras.igGrad_tptr;
    const label* lPtr                   = paras.lPtr;
    const label* uPtr                   = paras.uPtr;
    const label fSize                   = paras.fSize;
    const label cSize                   = paras.cSize;
    const label vector_size             = paras.vector_size;
    const label vector_offset           = paras.vector_offset;
    const label tensor_size             = paras.tensor_size;
    const label tensor_offset           = paras.tensor_offset;
    const struct rowSubsection** secs   = paras.secs;
    const label secNumInseg             = paras.secNumInseg;
    const label colRoundNum             = paras.colRoundNum;

    const label cpeId = athread_get_id(-1);
    //__thread_local char data[48000];
    label secIndex;
    /* calculate the row projection first */
    for(secIndex=0; secIndex< secNumInseg; secIndex++)
    {
        const struct rowSubsection subsection = secs[cpeId][secIndex];
        const label faceStart  = subsection.faceStart;
        const label nFace      = subsection.nFaces;

        /* get the row index and calculate row start and number */
        char* dataEnd = (char*)data;

        const label* owner = (const label*)ALIGNED(dataEnd);
        dataEnd = owner + nFace;

        volatile unsigned long get_reply=0;
        dma_desc gv=0;
        A_DMA_GET_SET( gv, PE_MODE, nFace*sizeof(label), &get_reply);
        A_DMA_GET_RUN( gv, &lPtr[faceStart], owner);
        DMA_Wait(&get_reply, 1);

        const label rowStart = owner[0];
        const label rowNum   = owner[nFace-1] - owner[0] + 1;

        // ldm_malloc memory
        const char* Sf = (const char*)ALIGNED(dataEnd);
        dataEnd = Sf + nFace*vector_size;

        const char* issf = (const char*)ALIGNED(dataEnd);
        dataEnd = issf + nFace*vector_size;

        char* igGrad = (char*)ALIGNED(dataEnd);
        dataEnd = igGrad + rowNum*tensor_size;

        if((label)( dataEnd - data) >= 48000)
        {
            CG_signals[cpeId] = 4;
            // core group must syncronize for colRoundNum
            // times in case normal cores block after bad
            // cores return
            label round;
            for(round = 0; round < secNumInseg*colRoundNum; round++)
                athread_syn( ARRAY_SCOPE, 0xFFFF);
            return 0;
        }


        // get the face fields
        dma_desc gv0=0, gv1=0, gv2=0;
        volatile unsigned long gReply0=0, gReply1=0, gReply2=0;

        A_DMA_GET_SET( gv0, PE_MODE, nFace*vector_size,    &gReply0);
        A_DMA_GET_SET( gv1, PE_MODE, nFace*vector_size, &gReply1);
        A_DMA_GET_SET( gv2, PE_MODE, rowNum*tensor_size, &gReply2);
        A_DMA_GET_RUN( gv0, &Sf_vptr[faceStart*vector_size],    Sf);
        A_DMA_GET_RUN( gv1, &issf_vptr[faceStart*vector_size],   issf);
        A_DMA_GET_RUN( gv2, &igGrad_tptr[rowStart*tensor_size], igGrad);

        DMA_Wait(&gReply0, 1);
        DMA_Wait(&gReply1, 1);
        DMA_Wait(&gReply2, 1);

        label faceIndex;
        for(faceIndex=0; faceIndex<nFace; faceIndex++){
            scalar* vSf     = &Sf[ faceIndex*vector_size + vector_offset ];
            scalar* vissf   = &issf[ faceIndex*vector_size + vector_offset ];
            scalar* vigGrad = &igGrad[ (owner[faceIndex]-rowStart)*tensor_size + tensor_offset ];

            vigGrad[0] += vSf[0]*vissf[0];
            vigGrad[1] += vSf[0]*vissf[1];
            vigGrad[2] += vSf[0]*vissf[2];
            vigGrad[3] += vSf[1]*vissf[0];
            vigGrad[4] += vSf[1]*vissf[1];
            vigGrad[5] += vSf[1]*vissf[2];
            vigGrad[6] += vSf[2]*vissf[0];
            vigGrad[7] += vSf[2]*vissf[1];
            vigGrad[8] += vSf[2]*vissf[2];
        }


        // add the tmp igGrad back
        dma_desc pv=0;
        volatile unsigned long pReply=0;
        A_DMA_PUT_SET( pv, PE_MODE, rowNum*tensor_size, &pReply);
        A_DMA_PUT_RUN( pv, igGrad, &igGrad_tptr[rowStart*tensor_size]);
        DMA_Wait(&pReply, 1);

    }


    /* calculate the column projection */
    for(secIndex=0; secIndex<secNumInseg; secIndex++)
    {

        const struct rowSubsection subsection = secs[cpeId][secIndex];
        const label faceStart  = subsection.faceStart;
        const label nFace      = subsection.nFaces;
        const label colRound   = subsection.colRound;
        const label nSecs      = subsection.nSecs;
//        const label * colSAndC = subsection.colStartsAndCounts;
		const label colSAndC[nSecs*2];
		dma_desc gv=0;
		volatile unsigned long gReply=0;
		A_DMA_GET_SET(gv,PE_MODE,nSecs*sizeof(label)*2,&gReply);
		A_DMA_GET_RUN(gv,subsection.colStartsAndCounts,&colSAndC[0]);
		DMA_Wait(&gReply,1);


        label totalCols = 0;
        label colSecIndex = nSecs;
        while(colSecIndex != 0){
            totalCols += colSAndC[2*colSecIndex-1];
            colSecIndex--;
        }

        char* dataEnd = (char*)data;

        const char* Sf = (const char*)ALIGNED(dataEnd);
        dataEnd = Sf + nFace*vector_size;

        const char* issf = (const char*)ALIGNED(dataEnd);
        dataEnd = issf + nFace*vector_size;

        const label* neibs = (const label*)ALIGNED(dataEnd);
        dataEnd = neibs + nFace;

        dma_desc gv0=0, gv1=0, gv2=0;
        volatile unsigned long gReply0=0, gReply1=0, gReply2=0;

        A_DMA_GET_SET( gv0, PE_MODE, nFace*vector_size,   &gReply0);
        A_DMA_GET_SET( gv1, PE_MODE, nFace*vector_size,   &gReply1);
        A_DMA_GET_SET( gv2, PE_MODE, nFace*sizeof(label), &gReply2);
        A_DMA_GET_RUN( gv0, &Sf_vptr[faceStart*vector_size],  Sf);
        A_DMA_GET_RUN( gv1, &issf_vptr[faceStart*vector_size], issf);
        A_DMA_GET_RUN( gv2, &uPtr[faceStart],                 neibs);

        scalar* tmp_igGrad = (scalar*)ALIGNED(dataEnd);
        dataEnd = tmp_igGrad + totalCols*9;

        char* igGrad = (char*)ALIGNED(dataEnd);
        dataEnd = igGrad + 1024; //(unsigned) (0.4*nFace)*tensor_size;
        if((label)( dataEnd - data) >= 48000)
        {
            CG_signals[cpeId] = 5;
            // core group must syncronize for colRoundNum
            // times in case normal cores block after bad
            // cores return
            label round;
            for(
                    round = 0;
                    round < colRoundNum*(secNumInseg-secIndex);
                    round++
                )
                athread_syn( ARRAY_SCOPE, 0xFFFF);
            return 0;
        }
        memset(tmp_igGrad, 0, totalCols*9*sizeof(scalar));

        DMA_Wait(&gReply0, 1);
        DMA_Wait(&gReply1, 1);
        DMA_Wait(&gReply2, 1);

        label faceIndex;
        for(faceIndex=0; faceIndex<nFace; faceIndex++){
            const scalar * vSf   = &Sf[ (faceIndex)*vector_size + vector_offset ];
            const scalar * vissf = &issf[ (faceIndex)*vector_size + vector_offset ];

            const label col_global = neibs[faceIndex];
            label col_local = 0;
            label secI = 0;
            while( col_global < colSAndC[2*secI] ||
                   col_global >= colSAndC[2*secI+1]+colSAndC[2*secI]
                 )
            {
                if( secI < nSecs ){
                    col_local += colSAndC[2*secI+1];
                    secI++;
                }
                else{
                        CG_signals[cpeId] = 6;
                        // core group must syncronize for colRoundNum
                        // times in case normal cores block after bad
                        // cores return
                        label round;
                        for(
                                round = 0;
                                round < colRoundNum*(secNumInseg-secIndex);
                                round++
                            )
                            athread_syn( ARRAY_SCOPE, 0xFFFF);
                        return 0;
                }
            }

            col_local += col_global - colSAndC[2*secI];

            tmp_igGrad[col_local*9 + 0] -= vSf[0]*vissf[0];
            tmp_igGrad[col_local*9 + 1] -= vSf[0]*vissf[1];
            tmp_igGrad[col_local*9 + 2] -= vSf[0]*vissf[2];
            tmp_igGrad[col_local*9 + 3] -= vSf[1]*vissf[0];
            tmp_igGrad[col_local*9 + 4] -= vSf[1]*vissf[1];
            tmp_igGrad[col_local*9 + 5] -= vSf[1]*vissf[2];
            tmp_igGrad[col_local*9 + 6] -= vSf[2]*vissf[0];
            tmp_igGrad[col_local*9 + 7] -= vSf[2]*vissf[1];
            tmp_igGrad[col_local*9 + 8] -= vSf[2]*vissf[2];
        }

        // communicate and write back according colRound
        label round;
        for(round = 0; round < colRoundNum; round++){
            athread_syn( ARRAY_SCOPE, 0xFFFF);
            if( round == colRound )
            {
                label secI;
                label tmpStart = 0;
                for( secI = 0; secI < nSecs; secI++){
                    label colStart = colSAndC[2*secI];
                    label colNum   = colSAndC[2*secI+1];

                    dma_desc gv3=0;
                    volatile unsigned long gReply3=0;
                    A_DMA_GET_SET( gv3, PE_MODE, colNum*tensor_size, &gReply3);
                    A_DMA_GET_RUN( gv3, &igGrad_tptr[colStart*tensor_size], igGrad);
                    DMA_Wait(&gReply3, 1);

                    label cellIndex;
                    for(cellIndex=0; cellIndex<colNum; cellIndex++)
                    {
                        scalar* vc = &igGrad[cellIndex*tensor_size + tensor_offset];
                        vc[0] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 0];
                        vc[1] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 1];
                        vc[2] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 2];
                        vc[3] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 3];
                        vc[4] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 4];
                        vc[5] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 5];
                        vc[6] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 6];
                        vc[7] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 7];
                        vc[8] += tmp_igGrad[ (cellIndex+tmpStart)*9 + 8];
                    }

                    dma_desc pv=0;
                    volatile unsigned long pReply=0;
                    A_DMA_PUT_SET( pv, PE_MODE, colNum*tensor_size, &pReply);
                    A_DMA_PUT_RUN( pv, igGrad, &igGrad_tptr[colStart*tensor_size]);
                    DMA_Wait(&pReply, 1);

                    tmpStart += colNum;
                }
            }
        }
    }
    return 0;
}

int swGgradfDivide_slave(swGgradfDivide_paras * paras)
{
	char data[48000];

    /* Global parameters */
    swGgradfDivide_paras para = *paras;
    scalar* igGrad = para.igGrad;
    const scalar* volume = para.volume;
    const label NCell    = para.NCell;
    const label size     = para.size;

    volatile label cpeID = athread_get_id(-1);

    /* decompsing parameters */
    volatile label segStart = ((NCell+63)/64)*cpeID;
    // take care of the array tail
    volatile label segSpan  =
        ((NCell+63)/64 < NCell-segStart) ?
        ((NCell+63)/64) : (NCell-segStart);
    volatile label segEnd   = segStart + segSpan;
    volatile label secSpan = (48000-2*32)/sizeof(scalar)/(size+1);

    /* Local parameters */
    //__thread_local char data[48000];
    char* dataEnd = data;
    scalar* igGrad_local = (scalar*)(ALIGNED(dataEnd));
    dataEnd = (char*)igGrad_local + secSpan*sizeof(scalar)*size;
    scalar* volume_local = (scalar*)(ALIGNED(dataEnd));
    dataEnd = (char*)volume_local + secSpan*sizeof(scalar);

    /* check whether the local Data exceeds bound */
    if( (label)(dataEnd - data) > 48000 )
    {
        CG_signals[cpeID] = 7;
        return 0;
    }

    volatile label secIter;
    volatile label localIter;
    volatile label cellIter;
    for(secIter = segStart; secIter<segEnd; secIter+=secSpan)
    {
        /* take care of segment tails */
        volatile label realSpan =
            (secSpan > segEnd-secIter) ? (segEnd-secIter) : (secSpan);

        /* get the local data */
        dma_desc gv1=0, gv2=0;
        volatile unsigned long get_reply1=0, get_reply2=0;
        A_DMA_GET_SET( gv1, PE_MODE, realSpan*sizeof(scalar)*size, &get_reply1);
        A_DMA_GET_SET( gv2, PE_MODE, realSpan*sizeof(scalar), &get_reply2);
        A_DMA_GET_RUN( gv1, &igGrad[secIter*size], igGrad_local);
        A_DMA_GET_RUN( gv2, &volume[secIter], volume_local);
        DMA_Wait(&get_reply1, 1);
        DMA_Wait(&get_reply2, 1);

        /* local core operations */
        for(localIter = 0; localIter<realSpan; localIter++ )
            for(cellIter = 0; cellIter<size; cellIter++ )
                igGrad_local[localIter*size+cellIter] /= volume_local[localIter];

        /* put back the local results */
        dma_desc pv=0;
        volatile unsigned long put_reply=0;
        A_DMA_PUT_SET( pv, PE_MODE, realSpan*sizeof(scalar)*size, &put_reply);
        A_DMA_PUT_RUN( pv, igGrad_local, &igGrad[secIter*size]);
        DMA_Wait(&put_reply, 1);
    }

    return 0;
}
