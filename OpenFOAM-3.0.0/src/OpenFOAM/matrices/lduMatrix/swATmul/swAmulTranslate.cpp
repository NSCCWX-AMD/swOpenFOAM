#include "swAmulTranslate.hpp"
extern "C" {
	#include "swAmul_host.h"
    #include <limits.h>
    #include "util.h"
}
#include <assert.h>
#include <memory.h>
#define SW_TRANS
#define DEBUG

void init_origin_location(const LABEL *matrix_size,
                          LABEL const *origin_location_ptr,
                          const LABEL *uPtr,
                          const LABEL *lPtr,
                          const LABEL nFaces,
                          const LABEL offset_size,
                          const LABEL nCell_block);

void translate_refill_location_host(refilltion_ptr parameter,
                                    const LABEL *uPtr,
                                    const LABEL *lPtr);

void init_each_block_psi_begin_end(amul_translate_array_ptr  data_array);
void assign_tasks(amul_translate_array_ptr data_array,
                  const LABEL Cell_Round_up_Block_size,
                  const LABEL *row_of_block_data_num);



refilltion translate_matrix_matrix(amul_translate_array_ptr &data_array,
                                   amul_para_ptr matrix,
                                   LABEL size,
                                   const LABEL offset_size)
{

    //the main function of translate for amul
    data_array = (amul_translate_array*)malloc(sizeof(amul_translate_array) * CORE_SIZE);

    const LABEL *uPtr = matrix->uPtr;
    const LABEL *lPtr = matrix->lPtr;
    const SCALAR *upperPtr = matrix->upperPtr;
    const SCALAR *lowerPtr = matrix->lowerPtr;
    amul_translate_ptr data;
    LABEL i=0;
    LABEL index;
    LABEL count = 0;
    LABEL column,row, index_lowwer;

    LABEL *matrix_size = (LABEL*) malloc(sizeof(LABEL) * size * size);
    LABEL *csr_row_index = (LABEL*) malloc( sizeof(LABEL) * (size + 1));
    const LABEL nFaces = matrix->nFaces;

    for(i = 0; i < CORE_SIZE; i++)
    {
        data_array[i].csr_row_size = size;
        data_array[i].csr_row_index = csr_row_index;
        data_array[i].nCells = matrix->nCells;
        data_array[i].offset_size = offset_size;
    }

    for(i = 0; i < size * size; i++)
    {
        matrix_size[i] = 0;
    }

    for(i = 0; i < size + 1; i++)
    {
        csr_row_index[i] = 0;
    }

    csr_row_index++;
    /////  count the number of data in each block and how much block which has data
    for(i = 0; i < matrix->nFaces; i++)
    {
        row = lPtr[i]>>offset_size;
        column = uPtr[i]>>offset_size;

        ////upper
        index = row * size + column;
        index_lowwer = row + column * size;
        matrix_size[index]++;
        count = matrix_size[index]==1 ? count+1:count;
        matrix_size[index_lowwer]++;
        count = matrix_size[index_lowwer]==1 ? count+1:count;
    }

    #define SIZE_TMP (sizeof(SCALAR)+(sizeof(LABEL)<<1))


    data = (amul_translate_ptr) malloc(sizeof(amul_translate) * count);
    for(i = 0; i < CORE_SIZE; i++)
    {
        data_array[i].data = data;
        data_array[i].size = count;
    }


    //////init origin location
    LABEL *origin_location_buffer = (LABEL*)malloc((sizeof(LABEL) * matrix->nFaces)) ;
    init_origin_location(matrix_size,
                         origin_location_buffer,
                         uPtr,
                         lPtr,
                         nFaces,
                         offset_size,
                         size);
    //////////////////////////


    //data_array->buffer = buffer;
    LABEL tmp, j;
    ////malloc the storage to store the translate matrix

    LABEL row_tmp,column_tmp;

    LABEL *upper_block_count = (LABEL*)malloc(sizeof(LABEL)*(size+1));

    for(i = 0; i < size +1; i++)
    {
        upper_block_count[i] = 0;
    }
    upper_block_count++;

    /////malloc memory for data
    for(i=0, j=0, row_tmp = 0, column_tmp = 0; i < size * size; i++)
    {

        tmp =  matrix_size[i];

        if(tmp > 0)
        {
            data[j].psi_begin = (column_tmp<<offset_size);
            data[j].Apsi_begin = (row_tmp<<offset_size);
            data[j].size = tmp;
            if(column_tmp == row_tmp)
            {
                char *buffer = (char*)malloc(sizeof(double) * tmp + sizeof(LABEL) * tmp);
                data[j].data = (SCALAR*)(buffer);

                data[j].uPtr = (LABEL*)(data[j].data + (tmp  ));
                data[j].lPtr = (LABEL*)(data[j].uPtr + (tmp >> 1 ));
            }
            else
            {
                char * buffer = (char*) malloc(sizeof(double) * tmp + sizeof(LABEL) * tmp * 2);
                data[j].data = (SCALAR*)(buffer);
                data[j].uPtr = (LABEL*)(data[j].data + (tmp ));
                data[j].lPtr = (LABEL*)(data[j].uPtr + (tmp ));
            }
            csr_row_index[row_tmp]++;
            //count the number of block in each row
            if(row_tmp < column_tmp)
            {
                upper_block_count[row_tmp]++;
            }
            j++;
        }
        column_tmp++;
        if(column_tmp==size)
        {
            column_tmp = 0;
            row_tmp++;
        }
    }

    for(i = 0; i < size - 1; i++)
    {
        upper_block_count[i + 1] += upper_block_count[i];
    }

    upper_block_count--;
    for(i = 0 ; i < size-1 ; i++ )
    {
        csr_row_index[i+1] += csr_row_index[i];

    }

//////////////init the number of data in each block row
    LABEL *upper_count = (LABEL*)malloc(sizeof(LABEL)*(size+1));

    get_block_row_end_index(lPtr, upper_count, nFaces, size, offset_size);
    //create lower ,  upper , and diag vector

////////////////////////////////////////////////
/////////////create lower blocks , upper blocks and diagnal blocks
    LABEL upper_block_num = upper_block_count[size];

    amul_translate_ptr upperPtr_block = (amul_translate_ptr) malloc(sizeof(amul_translate) * (upper_block_num * 2 + size));

    amul_translate_ptr lowerPtr_block = upperPtr_block + upper_block_num;
    amul_translate_ptr diagPtr_block = lowerPtr_block + upper_block_num;

    LABEL *lower_index = (LABEL*)malloc(sizeof(LABEL)*size);
    for(i = 0; i < size; i++)
    {
        lower_index[i] = 0 ;
    }
    for( i = 0 ; i < size ; i++)
    {
//      diagPtr_block[i].size = -1;
        diagPtr_block[i].size = 0 ;
    }

    LABEL w,w2;
    j = 0;
    amul_translate_ptr amul_translate_ptr_tmp;

    LABEL *lower_block_count = upper_block_count;
    for(i = 0; i < count; i++)
    {
        amul_translate_ptr_tmp = &(data[i]);
        if( amul_translate_ptr_tmp->psi_begin  > (amul_translate_ptr_tmp->Apsi_begin))
        {
            upperPtr_block[j++] = data[i];
        }
        else
        {
            if( amul_translate_ptr_tmp->psi_begin == amul_translate_ptr_tmp->Apsi_begin)
            {
                diagPtr_block[(amul_translate_ptr_tmp->Apsi_begin)>>offset_size] = data[i];
            }
            else
            {
                w =  amul_translate_ptr_tmp->psi_begin >>offset_size;
                w2 = lower_block_count[w] + lower_index[w];
                lowerPtr_block[ w2 ] = data[i];
                lower_index[w]++;
            }
        }
    }
//////////////////////////

///////////////////////////////////////init refill

    refilltion refill;
    refill.lowerPtr_block = lowerPtr_block;
    refill.upperPtr_block = upperPtr_block;
    refill.diagPtr_block = diagPtr_block;
    refill.upper_count = upper_count;
    refill.upper_block_count = upper_block_count;
    refill.lowerPtr = (double*)lowerPtr;
    refill.upperPtr = (double*)upperPtr;
    refill.nCell_block = size;
    refill.origin_location = origin_location_buffer;
    ///////fill data to memory

    translate_refill_location_host(&refill, uPtr, lPtr);

    // translate_refill_data( &refill );
    translate_refill_data_host(&refill);

    init_each_block_psi_begin_end(data_array);

    ///init row_begin end ,assign task to each duputy core
    assign_tasks(data_array,size,csr_row_index-1);
    free(matrix_size);
    free(lower_index);

    return refill;

    #undef SIZE_TMP
}

void assign_tasks(amul_translate_array_ptr data_array,
                  const LABEL Cell_Round_up_Block_size,
                  const LABEL *row_of_block_data_num)
{
    const LABEL each_core_has_row = ROUNDING_UP(Cell_Round_up_Block_size, DEPUTY_CORE_OFFSET);

    for(register LABEL i = 0; i < CORE_SIZE; i++)
    {
        const LABEL i_muil_each_core_has_row = i * each_core_has_row;
        data_array[i].block_row_begin = i_muil_each_core_has_row < Cell_Round_up_Block_size ? i_muil_each_core_has_row : Cell_Round_up_Block_size;
        data_array[i].block_row_end = i_muil_each_core_has_row + each_core_has_row < Cell_Round_up_Block_size ? i_muil_each_core_has_row + each_core_has_row : Cell_Round_up_Block_size;

        data_array[i].data = data_array[i].data + row_of_block_data_num[data_array[i].block_row_begin];

        data_array[i].block_num = row_of_block_data_num[data_array[i].block_row_end] - row_of_block_data_num[data_array[i].block_row_begin];
    }
}

#define ISDIAGNAL(block_tmp_ptr) ((block_tmp_ptr)->Apsi_begin == ((block_tmp_ptr)->psi_begin))

void init_each_block_psi_begin_end(amul_translate_array_ptr data_array)
{

    //this function is used to init the psi_begin psi_end in each block
    LABEL block_num = data_array->size;
    amul_translate_ptr blocks_ptr = data_array->data;

    for(LABEL i = 0; i < block_num; i++)
    {
        amul_translate_ptr block_tmp = blocks_ptr + i;

        const LABEL *uPtr = block_tmp->uPtr;
        register LABEL psi_begin = uPtr[0];
        register LABEL psi_end = uPtr[0];
        if(ISDIAGNAL(block_tmp))
        {
            const LABEL *lPtr = block_tmp->lPtr;
            const LABEL block_size = block_tmp->size >>1;

            for(LABEL j = 1; j < block_size; j++)
            {
                psi_begin = psi_begin > uPtr[j] ? uPtr[j] : psi_begin;
                psi_end = psi_end < uPtr[j] ? uPtr[j] : psi_end;
            }
            for(LABEL j = 0; j < block_size; j++)
            {
                psi_begin = psi_begin > lPtr[j] ? lPtr[j] : psi_begin;
                psi_end = psi_end < lPtr[j] ? lPtr[j] : psi_end;
            }
        }
        else
        {
            const LABEL block_size = block_tmp->size;

            for(LABEL j = 1; j < block_size; j++)
            {
                psi_begin = psi_begin > uPtr[j] ? uPtr[j] : psi_begin;
                psi_end = psi_end < uPtr[j] ? uPtr[j] : psi_end;
            }
        }

        block_tmp->psi_begin = psi_begin;
        block_tmp->psi_end = psi_end;
    }
}
#undef ISDIAGNAL



void translate_refill_data_host(refilltion_ptr parameter )
{
    //this function is used to init data(upperPtr,lowerPtr) by using the table(origin location)
    //conrefilltion refill = * parameter;
    amul_translate_ptr lowwerPtr_block  = parameter->lowerPtr_block;
    amul_translate_ptr upperPtr_block = parameter->upperPtr_block;
    amul_translate_ptr diagPtr_block = parameter->diagPtr_block;
    const SCALAR * upperPtr = parameter->upperPtr;
    const SCALAR * lowwerPtr = parameter->lowerPtr;

    const LABEL *origin_location = parameter->origin_location;

    LABEL *upper_block_count = parameter->upper_block_count;

    LABEL nCell_block = parameter->nCell_block;

    LABEL location_tmp;
    const LABEL * origin_location_tmp = origin_location;

    for(LABEL i = 0; i < nCell_block; i++)
    {
        //the order of init is 1st diagnal block , 1st row block , 1st column block
        //and then 2nd ...

        const LABEL row_block_end = upper_block_count[i + 1];

        const LABEL row_block_begin = upper_block_count[i];

        const LABEL calculate_size = diagPtr_block[i].size >> 1;
        //diagnal
        for (LABEL j = 0; j < calculate_size; ++j)
        {
            const LABEL location_tmp = origin_location_tmp[j];

            diagPtr_block[i].data[j] = upperPtr[location_tmp];
            diagPtr_block[i].data[j+calculate_size] = lowwerPtr[location_tmp];
        }
        origin_location_tmp += calculate_size;

        for(LABEL j = row_block_begin; j < row_block_end; j++)
        {
            ///upper
            const LABEL calculate_size = upperPtr_block[j].size;

            location_tmp = origin_location_tmp[0];
            for(LABEL k = 0; k < calculate_size; k++)
            {
                location_tmp = origin_location_tmp[k];

                upperPtr_block[j].data[k] = upperPtr[location_tmp];
            }
            //lowwer
            location_tmp = origin_location_tmp[0];
            for(LABEL k = 0; k < calculate_size; k++)
            {
                location_tmp = origin_location_tmp[k];

                lowwerPtr_block[j].data[k] = lowwerPtr[location_tmp];
            }
            origin_location_tmp += calculate_size;
        }
    }
}



void translate_refill_location_host(refilltion_ptr parameter, const LABEL *uPtr, const LABEL *lPtr)
{

    //this function is used to init uPtr and lPtr by using the table(origin location)
    amul_translate_ptr lowerPtr_block  = parameter->lowerPtr_block;
    amul_translate_ptr upperPtr_block = parameter->upperPtr_block;
    amul_translate_ptr diagPtr_block = parameter->diagPtr_block;

    const LABEL *origin_location = parameter->origin_location;

    LABEL *upper_block_count = parameter->upper_block_count;

    LABEL nCell_block = parameter->nCell_block;



    LABEL location_tmp;
    const LABEL *origin_location_tmp = origin_location;

    for(LABEL i = 0; i < nCell_block; i++)
    {
        //the order of init is 1st diagnal block , 1st row block , 1st column block
        //and then 2nd ...

        const LABEL row_block_end = upper_block_count[i + 1];

        const LABEL row_block_begin = upper_block_count[i];

        const LABEL calculate_size = diagPtr_block[i].size >> 1;
        //diagnal
        for (LABEL j = 0; j < calculate_size; ++j)
        {
            const LABEL location_tmp = origin_location_tmp[j];

            diagPtr_block[i].lPtr[j] = lPtr[location_tmp];
            diagPtr_block[i].uPtr[j] = uPtr[location_tmp];
        }
        origin_location_tmp += calculate_size;

        for(LABEL j = row_block_begin; j < row_block_end; j++)
        {
            ///upper
            const LABEL calculate_size = upperPtr_block[j].size;

            location_tmp = origin_location_tmp[0];
            for(LABEL k = 0; k < calculate_size; k++)
            {
                location_tmp = origin_location_tmp[k];
                upperPtr_block[j].lPtr[k] = lPtr[location_tmp];
                upperPtr_block[j].uPtr[k] = uPtr[location_tmp];
            }
            //lowwer
            location_tmp = origin_location_tmp[0];
            for(LABEL k = 0; k < calculate_size; k++)
            {
                location_tmp = origin_location_tmp[k];
                lowerPtr_block[j].lPtr[k] = uPtr[location_tmp];
                lowerPtr_block[j].uPtr[k] = lPtr[location_tmp];
            }
            origin_location_tmp += calculate_size;
        }
    }
}


void init_origin_location(const LABEL *matrix_size,
                          LABEL const *origin_location_ptr,
                          const LABEL *uPtr,
                          const LABEL *lPtr,
                          const LABEL  nFaces,
                          const LABEL  offset_size,
                          const LABEL  nCell_block)
{

    //init a table from origin data to data translated
    const LABEL nCell_block_to_1_sum = (nCell_block * (nCell_block + 1)) >>1;
	LABEL **matrix_origin_buffer_ptr_ptr = (LABEL**)malloc(sizeof(LABEL*)*nCell_block_to_1_sum);
    LABEL *matrix_origin_buffer_ptr_ptr_index = (LABEL*)malloc(sizeof(LABEL)*nCell_block_to_1_sum);

    LABEL *origin_location_ptr_tmp = (LABEL*)origin_location_ptr;
    LABEL index_tmp = 0;
    for(LABEL i = 0; i < nCell_block_to_1_sum; i++)
    {
        matrix_origin_buffer_ptr_ptr_index[i] = 0;
    }

    for(LABEL i = 0; i < nCell_block; i++)
    {
        //for diagnal
        const LABEL row_base = i * nCell_block;
        matrix_origin_buffer_ptr_ptr[index_tmp] = origin_location_ptr_tmp;
        origin_location_ptr_tmp += (matrix_size[ row_base + i] >> 1);
        index_tmp++;
        //for upper
        for(LABEL j = i + 1; j < nCell_block; ++j,index_tmp++)
        {
            matrix_origin_buffer_ptr_ptr[index_tmp] = origin_location_ptr_tmp;
            origin_location_ptr_tmp += matrix_size[row_base + j];
        }
    }

    for(LABEL i = 0; i < nFaces; i++)
    {
        const LABEL uPtri = uPtr[i]>>offset_size;
        const LABEL lPtri = lPtr[i]>>offset_size;

        const LABEL origin_location_block_index = (((nCell_block + nCell_block + 1 - lPtri) * (lPtri )) >> 1) + uPtri - lPtri;

        const LABEL matrix_origin_buffer_ptr_ptr_offset = matrix_origin_buffer_ptr_ptr_index[origin_location_block_index];
        matrix_origin_buffer_ptr_ptr[origin_location_block_index][matrix_origin_buffer_ptr_ptr_offset] = i ;
        matrix_origin_buffer_ptr_ptr_index[origin_location_block_index]++;
    }

    free(matrix_origin_buffer_ptr_ptr_index);
    free(matrix_origin_buffer_ptr_ptr);
}


amul_translate_array_ptr translate_matrix_matrix_T(const amul_translate_array_ptr data_array)
{
    //change the order of original block vector. to make the block be arranged in ascending
    //order if Apsi begin in Tmul calculate process
    //input is the block vector used in Amul
    amul_translate_array_ptr data_array_T;
    data_array_T = (amul_translate_array*)malloc(sizeof(amul_translate_array) * CORE_SIZE);

    const LABEL block_num = data_array[0].size;
    const LABEL matrix_size_div_block  = data_array[0].csr_row_size;
    amul_translate_ptr data_T = (amul_translate_ptr) malloc(sizeof(amul_translate) * block_num);
    const amul_translate_ptr  data = data_array[0].data ;

    for (LABEL i = 0; i < CORE_SIZE; ++i)
    {
        data_array_T[i] = data_array[i];

        data_array_T[i].data = data_T + (data_array[i].data - data);
    }

    const LABEL offset_size = data_array[0].offset_size;
    ///init data

    //init the number of block in each block .

    LABEL * each_row_block_num = (LABEL*)malloc(sizeof(LABEL) * matrix_size_div_block + sizeof(LABEL));
    memset(each_row_block_num, 0 , sizeof(LABEL) * matrix_size_div_block  + sizeof(LABEL));

    for (LABEL i = 0; i < data_array_T[0].size; ++i)
    {
        each_row_block_num[ (data[i].psi_begin>>offset_size) + 1 ]++;
    }
    ///add the former block number
    for (LABEL i = 1; i <= data_array_T[0].csr_row_size; ++i)
    {
        each_row_block_num[i] += each_row_block_num[i-1];
    }

    LABEL *each_row_block_index = (LABEL*)malloc(sizeof(LABEL) * data_array_T[0].csr_row_size);
    memset(each_row_block_index, 0, sizeof(LABEL) * (data_array_T[0].csr_row_size ));

    ///init block vector
    for (LABEL i = 0; i < block_num; ++i)
    {
        const LABEL block_row_index = data[i].psi_begin >> offset_size;
        const LABEL row_base = each_row_block_num[block_row_index];
        const LABEL index = row_base + each_row_block_index[block_row_index];
        data_T[ index ].size  = data[i].size;
        data_T[index].data = data[i].data;

        data_T[index].Apsi_begin = block_row_index << offset_size;
        data_T[index].psi_begin = data[i].Apsi_begin;
        data_T[index].uPtr = data[i].lPtr;
        data_T[index].lPtr = data[i].uPtr;

        each_row_block_index[block_row_index]++;
    }

    //init psi begin and end in block vector
    init_each_block_psi_begin_end(data_array_T);

    free(each_row_block_index);
    free(each_row_block_num);
    return data_array_T;
}
