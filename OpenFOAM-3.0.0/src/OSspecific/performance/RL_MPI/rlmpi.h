/* 
 * File:   prgrpcgc_slave.h
 * Author: Geng Chen (gengchn@gmail.com)
 * 
 * Created on December 17, 2017, 4:16 PM
 */
#ifndef FILTE_PRGPCGC_SLAVE_H
#define FILTE_PRGPCGC_SLAVE_H
#include <slave.h>
#include "RlmpiSharedType.h"
#define MaxNPackages 450
#define MaxNCycle 200
#define MaxNElm 35
#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)

#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst):"memory")
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst):"memory")
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var)::"memory")
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var)::"memory")

#define REG_SIMD_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst):"memory")
#define REG_SIMD_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst):"memory")
#define REG_SIMD_GETR(var) asm volatile ("getr %0\n":"=r"(var)::"memory")
#define REG_SIMD_GETC(var) asm volatile ("getc %0\n":"=r"(var)::"memory")
#define ROWSYN  athread_syn(ROW_SCOPE,0xff)
#define COLSYN  athread_syn(COL_SCOPE,0xff)
#define ALLSYN  athread_syn(ARRAY_SCOPE,0xffff)

#if USE_DYNAMIC_MEM_INDICE==1
extern __thread_local_fix volatile __attribute__ ((aligned(32))) Pack *_sPacks;
extern __thread_local_fix volatile __attribute__ ((aligned(32))) Pack *_rPacks;
extern __thread_local_fix Pack *_sPacks_same_col;
extern __thread_local_fix Pack *_sPacks_same_row;

extern __thread_local volatile __attribute__ ((aligned(32))) int8LDM* _putr_schedules;
extern __thread_local volatile __attribute__ ((aligned(32))) int8LDM* _getrputc_schedules;
extern __thread_local volatile __attribute__ ((aligned(32))) int8LDM* _getc_schedules;

extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _putr_schedules_same_row;
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _getr_schedules_same_row;

extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _putc_schedules_same_col;
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _getc_schedules_same_col;
#else
extern __thread_local_fix Pack _sPacks[MaxNPackages];
extern __thread_local_fix Pack _rPacks[MaxNPackages];
extern __thread_local_fix Pack *_sPacks_same_col;
extern __thread_local_fix Pack *_sPacks_same_row;
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putr_schedules[MaxNCycle];
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getrputc_schedules[MaxNCycle];
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getc_schedules[MaxNCycle];
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putr_schedules_same_row[MaxNCycle];
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getr_schedules_same_row[MaxNCycle];
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putc_schedules_same_col[MaxNCycle];
extern __thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getc_schedules_same_col[MaxNCycle];
#endif


extern __thread_local_fix volatile Table _table_ldm;
extern __thread_local_fix int _total_send_pcg;
extern __thread_local_fix int _total_recv_pcg;
extern __thread_local_fix int _nCycle;
extern __thread_local_fix int _nCycleSameCol;
extern __thread_local_fix volatile int _nCycleSameRow;

extern __thread_local int _get_reply, _put_reply;

extern inline void TransformPackage3(const Pack*sPacks, Pack*rPacks);
extern inline void TransformSameColumnPackage(const Pack*sPacks, Pack*rPacks);
extern inline void TransformSameRowPackage(const Pack*sPacks, Pack*rPacks);
extern inline void load_reg_mpi_init_data(Schedule * reg_data);
extern inline void REG_SIMD_GETR_PUTC();
#endif
