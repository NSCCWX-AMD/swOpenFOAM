/* 
 * File:   DataType.h
 * Author: Geng Chen (gengchn@gmail.com)
 * 
 * Created on December 17, 2017, 4:16 PM
 */
// This header define the data type used by both slave and host
#ifndef FILE_RLMPI_SHARED_TYPE_H
#define FILE_RLMPI_SHARED_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif

#define USE_DYNAMIC_MEM_INDICE 1
#define USE_DYNAMIC_MEM 1

#define MaxNPackages 450
#define MaxNCycle 200
#define MaxNElm 35
typedef int ThreadID;
#include <stdio.h>   
//#define FatalError(s) {                                             \
//  printf("Fatal error '%s' at %s:%d\n",s,__FILE__,__LINE__);        \
//  abort(); }


    typedef unsigned char int8LDM;
    typedef float sReal;
    typedef double dReal;
    typedef short int16LDM;

    typedef struct {
        int nPUTR;
        int nGETC;
        int nGETR_PUTC;

        int nGetrSameRow; //for same row communication
        int nPutrSameRow; //for same row communication

        int nGetcSameCol; //for same Col communication
        int nPutcSameCol; //for same Col communication

    } __attribute__ ((aligned(32))) Table; // 

    typedef struct {
        sReal data[6]; //data 32*6=192 bit
        unsigned res_pos : 16;
        unsigned src_id : 8;
        unsigned dst_id : 8;
        unsigned indM : 16; // convenient variable
        unsigned indP : 16; // convenient variable
    } __attribute__ ((aligned(32))) Pack; // 256 bit or 32 byte

    typedef struct {
        Table table[64];
        Pack *package[64];

        int nCycle;
        int8LDM *putr_schedules[64];
        int8LDM *getrputc_schedules[64];
        int8LDM *getc_schedules[64];

        int nCycleSameRow;
        int8LDM *putr_schedules_same_row[64];
        int8LDM *getr_schedules_same_row[64];

        int nCycleSameCol;
        int8LDM *putc_schedules_same_col[64];
        int8LDM *getc_schedules_same_col[64];
        int destroy;// If destroy all allocated memory
    } Schedule;

#define thread_mask 10000
#define mpi_mask 20000
#ifdef __cplusplus
}
#endif

#endif /* DATATYPE_H */

