//#include<stdio.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include "amulMacros.h"

#include <dma.h>

#define DIAG_BUFFER_SIZE 1024
#define UPPER_LOWER_BUFFER_SIZE (256)
#define UPPER_LOWER_BUFFER_SIZE_OFFSET (8)
#define MAX_BLOCK_NUM (1024)


void ATmul_Diag(amul_para_ptr amul_para_ptr_host)
{
    //this function is used to calculate the diag of matrix multiply vector
    // printf("amul used!!!\n");
    amul_para amul_slave;
    volatile int get_reply = 0, put_reply = 0;
    volatile int my_id;

    athread_get(PE_MODE,
                amul_para_ptr_host,
                &amul_slave,
                sizeof(amul_para),
                (int*) (&get_reply),
                0,0,0);

    while (get_reply != 1);

    SCALAR *amul_slave_ApsiPtr;
    SCALAR *amul_slave_diagPtr;
    my_id = athread_get_id(-1);

    LABEL size = ROUNDING_UP(amul_slave.nCells,DEPUTY_CORE_OFFSET);
    LABEL offset_vector = size * my_id;
    LABEL i,j;
    LABEL remaining = amul_slave.nCells - offset_vector;
    size = (remaining < size) ? remaining : size;

    if(size > 0)
    {
        SCALAR  amul_slave_ApsiPtr[DIAG_BUFFER_SIZE];
        SCALAR  amul_slave_diagPtr[DIAG_BUFFER_SIZE];
        LABEL remaining_buffer;
        for (i=0; i<size; i+=DIAG_BUFFER_SIZE, offset_vector+=DIAG_BUFFER_SIZE)
        {
            remaining_buffer = size - i;
            remaining_buffer = remaining_buffer < DIAG_BUFFER_SIZE ? remaining_buffer : DIAG_BUFFER_SIZE;
            //translate data
            get_reply = 0;
            athread_get(PE_MODE,
                        amul_slave.diagPtr + offset_vector,
                        amul_slave_diagPtr,
                        remaining_buffer * sizeof(SCALAR),
                        (int*)&get_reply,
                        0,0,0);
            athread_get(PE_MODE,
                        amul_slave.psiPtr + offset_vector,
                        amul_slave_ApsiPtr,
                        remaining_buffer * sizeof(SCALAR),
                        (int*)&get_reply,
                        0,0,0);
            //wait for data translate completed
            while(get_reply!=2);
            //begin calculate
            for(j=0; j<remaining_buffer; j++)
            {
                (amul_slave_ApsiPtr)[j] *= (amul_slave_diagPtr)[j];
            }
            put_reply = 0;
            athread_put(PE_MODE,
                        amul_slave_ApsiPtr,
                        amul_slave.ApsiPtr + offset_vector,
                        remaining_buffer * sizeof(SCALAR),
                        (int*)&put_reply,
                        0,0);
            while(put_reply!=1);
        }
    }
}

void ATmul_Upper_Lowwer(const amul_translate_array_ptr parameter_in)
{
    //this function is used for calculate the upper and lowwer data of matrix multiply vector
    //for both Amul and Tmul

    const int my_id = athread_get_id(-1);
    volatile int get_reply, put_reply=0;

    volatile amul_translate_array parameter;

    get_reply = 0;
    /*dma_desc dma_parameter;
    dma_set_op(&dma_parameter,DMA_GET);
    dma_set_reply(&dma_parameter,&get_reply);
    dma_set_stepsize(&dma_parameter,0);
    dma_set_mode(&dma_parameter,PE_MODE);
    dma_set_mask(&dma_parameter,0xff);
    dma_set_size(&dma_parameter,sizeof(amul_translate_array));

    dma(dma_parameter, (long)&(parameter_in[my_id]), (long)&parameter);

    dma_wait(&get_reply,1);     ///get the parameter*/
    athread_get(PE_MODE,
                &(parameter_in[my_id]),
                &parameter,
                sizeof(amul_translate_array),
                (int*) &get_reply,
                0,0,0);
    while(get_reply!=1);

    const LABEL block_boundary_size = (1<<parameter.offset_size);
    const LABEL offset_mask = (0xffffffff << parameter.offset_size);

    LABEL Apsi_size, psi_size, calculate_size, half_calculate_size;
    LABEL i, k;
    LABEL uPtr_slave[UPPER_LOWER_BUFFER_SIZE];
    LABEL lPtr_slave[UPPER_LOWER_BUFFER_SIZE];

    SCALAR ApsiPtr_slave[OFFSET_SIZE_TEN];
    SCALAR psiPtr_slave[OFFSET_SIZE_TEN];

    SCALAR data_slave[UPPER_LOWER_BUFFER_SIZE];  //the upper and lower data be store in here

    SCALAR *ApsiPtr_slave_shift; //used to shift the pointer of Apsi 'Apsi_begin' units
    SCALAR *psiPtr_slave_shift;  //similar to the above


    LABEL psi_begin_now, psi_end_now, Apsi_begin_now;

    amul_translate amul_translate_ptr_slave[MAX_BLOCK_NUM];    //calculate blocks

    const SCALAR *parameter_ApsiPtr = parameter.ApsiPtr;
    const SCALAR *parameter_psiPtr  = parameter.psiPtr;
    const LABEL parameter_offset_size = parameter.offset_size;
    const LABEL nCells = parameter.nCells;

    LABEL block_index;
    for (block_index=0; block_index<parameter.block_num; block_index+=MAX_BLOCK_NUM)
    {
        const LABEL block_num = parameter.block_num - block_index < MAX_BLOCK_NUM ? parameter.block_num - block_index : MAX_BLOCK_NUM ;

        get_reply = 0 ;
        athread_get(PE_MODE,
                    parameter.data + block_index,
                    amul_translate_ptr_slave,
                    block_num * sizeof(amul_translate),
                    (int*) &get_reply,
                    0,0,0);
        while(get_reply!=1);

        //init dma descriptor of psi
        /*dma_desc dma_psi;
        dma_set_op(&dma_psi,DMA_GET);
        dma_set_reply(&dma_psi,&get_reply);
        dma_set_stepsize(&dma_psi,0);
        dma_set_mode(&dma_psi,PE_MODE);
        dma_set_mask(&dma_psi,0xff);*/

        LABEL index;
        for (index=0; index<block_num;)
        {
            //for each line
            //get Apsi
            Apsi_begin_now = amul_translate_ptr_slave[index].Apsi_begin;

            ApsiPtr_slave_shift =  (SCALAR *)(ApsiPtr_slave - Apsi_begin_now);
            Apsi_size = (nCells - Apsi_begin_now);
            Apsi_size = Apsi_size<block_boundary_size? Apsi_size:block_boundary_size;

            get_reply =0;
            athread_get(PE_MODE,
                        (SCALAR*)(parameter_ApsiPtr + Apsi_begin_now),
                        (SCALAR*) ApsiPtr_slave,
                        Apsi_size *(sizeof(SCALAR)),
                        (int*) &get_reply,
                        0,0,0);
            while((get_reply)!=1); //has get Apsi


            const LABEL Apsi_begin_mask = Apsi_begin_now;

            for(; index<block_num&&(Apsi_begin_mask==(amul_translate_ptr_slave[index].Apsi_begin)); index++)
            {
                //for blocks in  line
                psi_begin_now = amul_translate_ptr_slave[index].psi_begin;
                psiPtr_slave_shift = (SCALAR *)(psiPtr_slave - psi_begin_now);
                psi_end_now = amul_translate_ptr_slave[index].psi_end;
                psi_size = (psi_end_now - psi_begin_now + 1);

                get_reply = 0;

                /*dma_set_size(&dma_psi, sizeof(SCALAR)*psi_size);
                dma(dma_psi, (long)(parameter_psiPtr + psi_begin_now), (long)psiPtr_slave);*/
                athread_get(PE_MODE,
                        parameter_psiPtr + psi_begin_now,
                        psiPtr_slave,
                        sizeof(SCALAR)*psi_size,
                        (int*) &get_reply,
                        0,0,0);
                while(get_reply!=1);

                calculate_size = amul_translate_ptr_slave[index].size;
                half_calculate_size = calculate_size >> 1;
                // while(get_reply!=1);

                //dma_wait(&get_reply,1);

                if((psi_begin_now & offset_mask) == Apsi_begin_now)
                {
                    //diagnal block
                    LABEL buffer_index;
                    for (buffer_index=0; buffer_index<half_calculate_size; buffer_index+=UPPER_LOWER_BUFFER_SIZE)
                    {
                        get_reply = 0 ;
                        const LABEL remaining_buffer = half_calculate_size - buffer_index < UPPER_LOWER_BUFFER_SIZE ? half_calculate_size - buffer_index : UPPER_LOWER_BUFFER_SIZE;

                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].data + half_calculate_size + buffer_index ) ,
                                    (void*)(data_slave),
                                    remaining_buffer * sizeof(SCALAR),
                                    (int*) &get_reply,
                                    0,0,0 );
                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].uPtr + buffer_index ),
                                    (void*)( uPtr_slave ),
                                    remaining_buffer * sizeof(LABEL),
                                    (int*) &get_reply,
                                    0,0,0 );
                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].lPtr + buffer_index ),
                                    (void*)(lPtr_slave),
                                    remaining_buffer * sizeof(LABEL),
                                    (int*) &get_reply,
                                    0,0,0 );
                        while(get_reply!=3);

                        //calcalate lower
                        for(k=0; k<remaining_buffer; k++)
                        {
                            ApsiPtr_slave_shift[uPtr_slave[k]] += data_slave[k] *
                                psiPtr_slave_shift[lPtr_slave[k]];
                        }

                        //////////
                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].data + buffer_index),
                                    (void*)(data_slave),
                                    remaining_buffer * sizeof(SCALAR),
                                    (int*) &get_reply,
                                    0,0,0);
                        while(get_reply!=4);
                        //has get all data needed
                        //calculate upper
                        for(k=0; k<remaining_buffer; k++)
                        {
                            ApsiPtr_slave_shift[lPtr_slave[k]] += data_slave[k] * psiPtr_slave_shift[uPtr_slave[k]];
                        }
                    }//endl of buffer
                }
                else
                {
                    int buffer_index;
                    for (buffer_index=0; buffer_index<calculate_size; buffer_index+=UPPER_LOWER_BUFFER_SIZE)
                    {
                        get_reply = 0;
                        const LABEL remaining_buffer = calculate_size - buffer_index < UPPER_LOWER_BUFFER_SIZE ? calculate_size - buffer_index : UPPER_LOWER_BUFFER_SIZE;

                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].data + buffer_index),
                                    (void*)(data_slave),
                                    remaining_buffer * sizeof(SCALAR),
                                    (int*) &get_reply,
                                    0,0,0);

                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].uPtr + buffer_index),
                                    (void*)(uPtr_slave),
                                    remaining_buffer * sizeof(LABEL),
                                    (int*) &get_reply,
                                    0,0,0);

                        athread_get(PE_MODE,
                                    (void*)(amul_translate_ptr_slave[index].lPtr + buffer_index),
                                    (void*)(lPtr_slave),
                                    remaining_buffer * sizeof(LABEL),
                                    (int*) &get_reply,
                                    0,0,0);
                        while(get_reply!=3);

                         //has get all data needed
                        for(k=0; k<remaining_buffer; k++)
                        {
                            ApsiPtr_slave_shift[lPtr_slave[k]] += data_slave[k] * psiPtr_slave_shift[uPtr_slave[k]];
                        }
                    }//end of buffer
                }//endl of else
            }//endl of row for

            //put Apsi to host
            put_reply = 0;
            athread_put(PE_MODE,
                        (SCALAR*) (ApsiPtr_slave),
                        (SCALAR*) (parameter_ApsiPtr + Apsi_begin_now),
                        Apsi_size * sizeof(SCALAR),
                        (int*) &put_reply,
                        0,0);   //put old data to host
            while(put_reply!=1);
        }//end of each block row
    }//end of block buffer cycle
}//end of fuction
