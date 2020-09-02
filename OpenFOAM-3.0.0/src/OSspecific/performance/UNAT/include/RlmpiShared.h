/* 
 * File:   SharedType.h
 * Author: Geng Chen (gengchn@gmail.com)
 * 
 * Created on December 17, 2017, 4:16 PM
 */
// This header define the data type used by both slave and host
#ifndef FILE_RLMPI_SHARED_H
#define FILE_RLMPI_SHARED_H

#define USE_DYNAMIC_MEM_INDICE 1
#define USE_DYNAMIC_MEM 1

#define MaxNPackages 1450
#define MaxNCycle 255

typedef int int32;
typedef unsigned char int8LDM;
typedef float sReal;
typedef double dReal;
typedef short int16LDM;
typedef int IndLDM;


#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>   
//#define FatalError(s) {                                             
//  printf("Fatal error '%s' at %s:%d\n",s,__FILE__,__LINE__);        
//  abort(); }

    typedef struct {
        int nPutrSkew;
        int nGetcSkew;
        int nGetrPutcSkew;

        int nGetrSameRow; //for same row communication
        int nPutrSameRow; //for same row communication

        int nGetcSameCol; //for same Col communication
        int nPutcSameCol; //for same Col communication

    } Table;

    typedef struct {
        sReal data[6]; //data 32*6=192 bit
        unsigned res_pos : 16;
        unsigned src_id : 8;
        unsigned dst_id : 8;
        unsigned cva : 16; // convenient variable a
        unsigned cvb : 16; // convenient variable b
    } __attribute__ ((aligned(32))) Pack; // 256 bit or 32 byte

    typedef struct {
        Table table[64];
        Pack *package[64];

        int8LDM* srcId_list[64];
        int8LDM* dstId_list[64];
        int16LDM* resPos_list[64];
        
        int16LDM* cva_list[64];
        int16LDM* cvb_list[64];

        int nCycle;
        int8LDM* putr_schedules[64];
        int8LDM* getrputc_schedules[64];
        int8LDM* getc_schedules[64];

        int nCycleSameRow;
        int8LDM* putr_schedules_same_row[64];
        int8LDM* getr_schedules_same_row[64];

        int nCycleSameCol;
        int8LDM *putc_schedules_same_col[64];
        int8LDM *getc_schedules_same_col[64];
    } RlmpiInfo;

void slave_initRlmpiInfo(RlmpiInfo* rlmpi_info);
void slave_destroyRlmpiInfo(RlmpiInfo* rlmpi_info);
#ifdef __cplusplus
}
#endif

#endif /* DATATYPE_H */

