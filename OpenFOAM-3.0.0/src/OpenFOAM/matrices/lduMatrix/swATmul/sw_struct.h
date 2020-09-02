#ifndef SW_STRUCT_H
#define SW_STRUCT_H

#define SCALAR double
#define LABEL int
#define DEPUTY_CORE_OFFSET (6)
//#define DEPUTY_CORE_OFFSET (1)
#define CORE_SIZE (1<<DEPUTY_CORE_OFFSET)
//#define CORE_SIZE 64

#define OFFSET_SIZE (8)
#define OFFSET_SIZE_TEN (256)


#define ALIGNED(addr) ((( ( (unsigned long)(addr) - 1)>>5)+1)<<5)


#define LESS_TARGET_OFFSET(data,target,offset) ((data>>offset<<offset)<target)
#define NUM_LEFT_OFF_NUM(a,b) ((a)<<(b))

#define ROUNDING_UP(divisor,offset) ((((divisor)-1)>>(offset))+1)
#define SET_OFFSET_ZERO(a,offset) ((a)>>offset<<offset)

#define ArraySize 57344
#define TRANSLATE_BUFFER  (1024)

#define SW_TMUL_COARSE_RELATION(matrix,w,x,interfaceBouCoeffs,interfaces,cmpt,amul_parameter_ptr,matrix_translate_ptr,coarseLevel)  \
    do  {\
            if((coarseLevel) < MAX_SW_USING_CORASE_LEVELS)\
            {\
                SW_Tmul(matrix, w, x, interfaceBouCoeffs, interfaces, cmpt, *(amul_parameter_ptr), (matrix_translate_ptr));\
            }\
            else\
            {\
                matrix.Tmul(w, x, interfaceBouCoeffs, interfaces, cmpt);\
            }\
        }while(0)

#define SW_AMUL_COARSE_RELATION_REDESIGNED_12(matrix,w,x,interfaceBouCoeffs,interfaces,cmpt,amul_parameter_ptr,matrix_translate_ptr,coarseLevel)  \
    do  {\
            if((coarseLevel) < MAX_SW_USING_CORASE_LEVELS)\
            {\
                SW_Amul(matrix, w, x, interfaceBouCoeffs, interfaces, cmpt, *(amul_parameter_ptr), (matrix_translate_ptr));\
            }\
            else\
            {\
                matrix.Amul(w, x, interfaceBouCoeffs, interfaces, cmpt);\
            }\
        }while(0)

#define SW_AMUL_COARSE_RELATION_REDESIGNED_11(matrix,w,x,interfaceBouCoeffs,interfaces,cmpt,amul_parameter_ptr,matrix_translate_ptr,coarseLevel)  \
    do  {\
                matrix.Amul(w, x, interfaceBouCoeffs, interfaces, cmpt);\
        }while(0)

#define SW_AMUL_COARSE_RELATION_REDESIGNED(matrix,w,x,interfaceBouCoeffs,interfaces,cmpt,amul_parameter_ptr,matrix_translate_ptr,coarseLevel)  \
    do  {\
            if(!coarseLevel)\
            {\
                swAmulRedesigned(matrix, w, x, interfaceBouCoeffs, interfaces, cmpt, *(amul_parameter_ptr), (matrix_translate_ptr));\
            }\
            else if((coarseLevel) < MAX_SW_USING_CORASE_LEVELS)\
            {\
                swGAMGAmulRedesigned(matrix, w, x, interfaceBouCoeffs, interfaces, cmpt, *(amul_parameter_ptr), (matrix_translate_ptr));\
            }\
            else\
            {\
                GAMGAmulRedesigned(matrix, w, x, interfaceBouCoeffs, interfaces, cmpt);\
            }\
        }while(0)


#define SW_ATMUL_COARSE_RELATION(matrix,w,wT,x,interfaceBouCoeffs,interfaceIntCoeffs,interfaces,cmpt,amul_parameter_ptr,matrix_translate_ptr,coarseLevel)  \
    do {                                                                                                                            \
           \
                SW_ATmul(matrix, w, wT, x, interfaceBouCoeffs, interfaceIntCoeffs, interfaces, cmpt,*(amul_parameter_ptr), (matrix_translate_ptr));        \
    }while(0)

typedef struct amul_translate_tmp {
    LABEL *psi_begin;
    LABEL *psi_end;
    //LABEL psi_end;
    LABEL Apsi_begin;
    LABEL size;
    SCALAR  *data;
    LABEL  *uPtr;
    LABEL  *lPtr;
    LABEL * origin_location;
} amul_translate_tmp,*amul_translate_tmp_ptr;


typedef struct amul_translate_refill {
    LABEL size;           //block size
    SCALAR  *data;      //data in block
} amul_translate_refill,*amul_translate_refill_ptr;

typedef struct amul_translate {
    LABEL psi_begin;
    LABEL psi_end;
    LABEL Apsi_begin;
    LABEL size;           //block size
    SCALAR  *data;      //data in block
    LABEL *uPtr;         //row index in block
    LABEL *lPtr;         //column index in block
} amul_translate,*amul_translate_ptr;
typedef struct amul_translate_array {
//    LABEL my_id ;
    LABEL size ;
    LABEL block_row_begin;
    LABEL block_row_end;
    LABEL nCells;
    LABEL csr_row_size;
    LABEL offset_size ;
    LABEL block_num;
    LABEL *csr_row_index;

    SCALAR * psiPtr;
    SCALAR * ApsiPtr;
    amul_translate_ptr data;
} amul_translate_array, *amul_translate_array_ptr;

typedef struct scaling_factor_para {
    LABEL size;
    SCALAR * source;
    SCALAR * field;
    SCALAR * Acf;
    SCALAR * D;
    SCALAR * scalingFactorNum ;
    SCALAR * scalingFactorDenom;
} scaling_factor_para, *scaling_factor_para_ptr;
typedef struct amul_para {
    SCALAR * ApsiPtr;
    SCALAR * psiPtr;
    SCALAR * diagPtr;
    SCALAR * lowerPtr;
    SCALAR * upperPtr;
    LABEL *lPtr;
    LABEL *uPtr;
    LABEL nFaces;
    LABEL nCells;
    LABEL *matrix_size;

} *amul_para_ptr,amul_para;
typedef struct refilltion {
    SCALAR * upperPtr;
    SCALAR * lowerPtr;
    //LABEL maxBlockSize;
//    LABEL * lPtr;
  //  LABEL * uPtr;
//    LABEL nFaces;
    LABEL nCell_block;
    LABEL * origin_location;
    //LABEL ** lowwer_origin_location;
    //LABEL ** upper_origin_location;
    //LABEL ** dia_origin_location ;
    LABEL * upper_block_count;
    LABEL * upper_count;
    amul_translate_ptr lowerPtr_block ;
    amul_translate_ptr upperPtr_block ;
    amul_translate_ptr diagPtr_block;
    //amul_translate_array_ptr data_array;
    //LABEL ** origin_location;
} refilltion , * refilltion_ptr;


void init_amul_para(amul_para_ptr parameter,SCALAR *diag,SCALAR* psi,SCALAR * Apsi,SCALAR * lowerPtr,SCALAR *upperPtr,LABEL *lPtr,LABEL *uPtr,LABEL NFACES,LABEL SIZE);

#endif
