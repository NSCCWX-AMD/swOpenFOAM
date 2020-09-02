
#ifndef REGSTER_COMMUNICATION_H
#define REGSTER_COMMUNICATION_H

#include <math.h>
#include "slave.h"
#include "simd.h"
#define ROW_SIZE 8
#define COLUMN_SIZE 8
#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define ROWSYN()  athread_syn(ROW_SCOPE,0xff)
#define COLSYN()  athread_syn(COL_SCOPE,0xff)
#define ALLSYN()  athread_syn(ARRAY_SCOPE,0xffff)


#endif

