#ifndef ONCHIPTRANS_H
#define ONCHIPTRANS_H

#include "swMacro.h" 
#include "slaveUtils.h"
#include "RlmpiShared.h"

typedef int oct_unit;
#ifndef ROW
#define ROW(id) ((id) >> 3)
#endif
#ifndef COL
#define COL(id) ((id) & 0x07)
#endif

#define PACKSIZE (24/sizeof(oct_unit))
typedef struct {
	oct_unit val[PACKSIZE]; //data 32*6=192 bit
	unsigned src     : 8;
	unsigned dest    : 8;
	unsigned valSize : 16; 
	unsigned indD    : 16;// convenient variable a
	unsigned indM    : 16;// convenient variable b
} __attribute__ ((aligned(32))) Pack_tmp; // 256 bit or 32 byte
typedef  Pack_tmp oct_pack;


#define UNPACK_DATA( pack, data )\
{\
	int ival;\
	for(ival = 0; ival < (pack).valSize; ival++)\
		(data)[ival] = (pack).val[ival];\
	data += (pack).valSize;\
}

//if(_MYID==ID && pack.src==10) printf("src: %d, valSize: %d, val: %.8lf\n",pack.src,pack.valSize,((double*) &pack.val[0])[0]); \

#define PACK_DATA( pack, data )\
{\
	int ival;\
	for(ival = 0; ival < (pack).valSize; ival++)\
		(pack).val[ival] = (data)[ival];\
	data += pack.valSize; \
}

#define ARRAY_PUTR( data, size, destId)\
{\
	swInt packSize = ((size) / PACKSIZE );\
	swInt restCount = \
		( (size) % PACKSIZE );\
	oct_pack pack;\
	pack.src = _MYID;\
	pack.dest = destId;\
	pack.valSize = PACKSIZE;\
	swInt col_dest = COL(destId); \
	swInt j;\
	for(j = 0; j < packSize; j++)\
	{\
		PACK_DATA(pack, data); \
		REG_PUTR( *(int256*) &pack , col_dest);\
	}\
	if(restCount != 0)\
	{\
		pack.valSize = restCount;\
		PACK_DATA(pack, data); \
		REG_PUTR( *(int256*) &pack , col_dest);\
	}\
}

#define ARRAY_PUTC( data, size, destId)\
{\
	swInt packSize = size / PACKSIZE;\
	swInt restCount = \
		( size % PACKSIZE );\
	oct_pack pack;\
	pack.src = _MYID;\
	pack.dest = destId;\
	pack.valSize = PACKSIZE;\
	swInt j;\
	for(j = 0; j < packSize; j++)\
	{\
		PACK_DATA( pack, data);\
		REG_PUTC( *(int256*) &pack, ROW(destId));\
	}\
	if(restCount != 0)\
	{\
		pack.valSize = restCount;\
		PACK_DATA( pack, data);\
		REG_PUTC( *(int256*) &pack, ROW(destId));\
	}\
}


// on chip transfer operator
#define OCT_MAX(target, src) \
{ (target) = (target) > (src) ? (target) : (src); }

#define OCT_MIN(target, src)\
{ (target) = (target) < (src) ? (target) : (src); }

#define OCT_PLUS(target, src) { (target) += (src); }

#define OCT_MULT(target, src) { (target) *= (src); }

// all reduce
#include "oct_allReduce.h"
//// reduce
//#include "oct_reduce.h"
// p2p transfer
#define p2pTrans(data, size, src, dest) \
{ \
}
// group transfer
//void groupTrans(oct_unit** sdata, // 64 to avoid localization
//				swInt snum, swInt* ssizes, swInt* dests, 
//				oct_unit** rdata, // 64 to avoid localization 
//				swInt rnum, swInt* rsizes, swInt* srcs)
//{
//}
#define groupTrans(sdata, ssizes, rdata, rsizes) \
{ \
	INIT_LDM_SPACE(2048); \
	swInt i; \
	int icol,irow,idx; \
	oct_unit **sendData,**recvData; \
	LDM_NEW(sendData,oct_unit*,BLOCKNUM64K); \
	LDM_NEW(recvData,oct_unit*,BLOCKNUM64K); \
	for(i=0;i<BLOCKNUM64K;i++) \
	{ \
		sendData[i] = (oct_unit*)sdata[i]; \
		recvData[i] = (oct_unit*)rdata[i]; \
	} \
	oct_pack pack;  \
	swInt rowID = ROW(_MYID); \
	swInt colID = COL(_MYID); \
	/****right down*****/ \
	swInt num_Get = 0; \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( ROW(i) < rowID \
					&& COL(i) < colID) \
			num_Get += (rsizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
	} \
	swInt num_switch[8] = {0}; \
	for(i = 0; i<BLOCKNUM64K; i++)  \
	{ \
		if( COL(i) > colID && \
					ROW(i) >= rowID ) \
		{ \
			num_switch[COL(i)] +=  \
				(ssizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
		} \
	} \
	oct_rowAllReduce(num_switch, swInt, 8, OCT_PLUS); \
	for(i = 0; i < num_Get ; i++) \
	{ \
		REG_GETC( *(int256*) &pack ); \
		UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(i = 0; i < num_switch[colID]; i++)  \
	{ \
		REG_GETR( *(int256*) &pack ); \
		if( pack.dest != _MYID ) \
			REG_PUTC( *(int256*) &pack, ROW(pack.dest)); \
		else \
			UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(icol=colID+1;icol<8;icol++) \
	{ \
		for(irow=rowID;irow<8;irow++) \
		{ \
			idx = icol+irow*8; \
			ARRAY_PUTR(sendData[idx], ssizes[idx], idx); \
		} \
	} \
	ALLSYN; \
	/****left down*****/ \
	num_Get = 0;  \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( ROW(i) < rowID \
					&& COL(i) > colID) \
			num_Get += (rsizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
	} \
	for(i=0;i<8;i++) { num_switch[i] = 0; } \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( COL(i) <= colID && \
					ROW(i) > rowID ) \
		{ \
			num_switch[ROW(i)] +=  \
				(ssizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
		} \
	} \
	oct_colAllReduce(num_switch, swInt, 8, OCT_PLUS); \
	for(i = 0; i < num_Get ; i++) \
	{ \
		REG_GETR( *(int256*) &pack ); \
		UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(i = 0; i < num_switch[rowID]; i++) \
	{ \
		REG_GETC( *(int256*) &pack ); \
		if( pack.dest != _MYID ) \
			REG_PUTR( *(int256*) &pack, COL(pack.dest)); \
		else \
			UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(irow=rowID+1;irow<8;irow++) \
	{ \
		for(icol=colID;icol>=0;icol--) \
		{ \
			idx = icol+irow*8; \
			ARRAY_PUTC(sendData[idx], ssizes[idx], idx); \
		} \
	} \
	ALLSYN; \
	/****right up*****/ \
	num_Get = 0; \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( ROW(i) > rowID \
					&& COL(i) < colID) \
			num_Get += (rsizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
	} \
	for(i=0;i<8;i++) { num_switch[i] = 0; } \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( COL(i) >= colID && \
					ROW(i) < rowID ) \
		{ \
			num_switch[ROW(i)] +=  \
				(ssizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
		} \
	} \
	oct_colAllReduce(num_switch, swInt, 8, OCT_PLUS);\
	for(i = 0; i < num_Get ; i++) \
	{ \
		REG_GETR( *(int256*) &pack ); \
		UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(i = 0; i < num_switch[rowID]; i++) \
	{ \
		REG_GETC( *(int256*) &pack ); \
		if( pack.dest != _MYID ) \
			REG_PUTR( *(int256*) &pack, COL(pack.dest)); \
		else \
			UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(irow=rowID-1;irow>=0;irow--) \
	{ \
		for(icol=colID;icol<8;icol++) \
		{ \
			idx = icol+irow*8; \
			ARRAY_PUTC(sendData[idx], ssizes[idx], idx); \
		} \
	} \
	ALLSYN; \
	/****left up*****/ \
	num_Get = 0; \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( ROW(i) > rowID \
					&& COL(i) > colID) \
			num_Get += (rsizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
	} \
	for(i=0;i<8;i++) { num_switch[i] = 0; } \
	for(i = 0; i<BLOCKNUM64K; i++) \
	{ \
		if( COL(i) < colID && \
					ROW(i) <= rowID ) \
		{ \
			num_switch[COL(i)] +=  \
				(ssizes[i] + PACKSIZE - 1) \
				/ PACKSIZE; \
		} \
	} \
	oct_rowAllReduce(num_switch, swInt, 8, OCT_PLUS);\
	for(i = 0; i < num_Get ; i++) \
	{ \
		REG_GETC( *(int256*) &pack ); \
		UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(i = 0; i < num_switch[colID]; i++) \
	{ \
		REG_GETR( *(int256*) &pack ); \
		if( pack.dest != _MYID ) \
			REG_PUTC( *(int256*) &pack, ROW(pack.dest)); \
		else \
			UNPACK_DATA(pack, recvData[pack.src]); \
	} \
	for(icol=colID-1;icol>=0;icol--) \
	{ \
		for(irow=rowID;irow>=0;irow--) \
		{ \
			idx = icol+irow*8; \
			ARRAY_PUTR(sendData[idx], ssizes[idx], idx); \
		} \
	} \
	ALLSYN; \
}

#endif //ONCHIPTRANS_H

















