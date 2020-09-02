#ifndef SMVM_H
#define SMVM_H
#include <stdlib.h>
#include <athread.h>
#include "sw_struct.h"
//#include "SMVM_routine.h"
//#define NONZERONUM 1516800
//#define CELLNUM 512000
#define BLOCKNUM64K 64
#define EPSILON 1.0e-4
#define DEBUG 0
#define MONITORIDX 21211
#define TEST 1
//extern "C" {
//	void SMVM_SP(){
//		athread_init();
//		athread_spawn(func,0);
//		athread_join();
//		athread_halt();
//	};
//};
//void SMVM_SP();
extern SLAVE_FUN(func)();

//struct LDUMatrix{
//	SCALAR* diag;
//	SCALAR* lower;
//	SCALAR* upper;
//	int* rowAddr;
//	int* colAddr;
//	int numCell;
//	int numEdge;
//};
//void SMVM_MP_BlockOrder(struct LDUMatrix Matrix, SCALAR *x, SCALAR *b_MLBSlaveReordered);
void SMVM_SP_BlockOrder(struct LDUMatrix Matrix);
void BlockOrderingMLB(struct LDUMatrix Matrix);
void initOwnNeiSendList();

//int max(int v1,int v2){return v1>v2? v1:v2;};
SCALAR *upperReordered,*lowerReordered,*diagReordered,*xReordered,*b_MLBSlave,*b_MLB,*b_MLBSlaveMP;
int *ownReordered,*neiReordered,*blockStarts,*cellStarts,*postCellOrder,*postEdgeOrder,*ownNeiSendList;
int cpeBlockNum,maxEdges,mtxBlockNum,mshBlockNum,maxCells,spIndex,maxBlock;

#endif
