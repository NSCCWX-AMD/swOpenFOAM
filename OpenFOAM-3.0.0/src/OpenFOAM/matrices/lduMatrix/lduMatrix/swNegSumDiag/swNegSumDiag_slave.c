#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include <assert.h>
#include <dma.h>
#include <unistd.h>
#include <simd.h>
#include <math.h>
#include "rlmpi.h"
#include "sw_struct.h"

#define BLOCKNUM64K 64
#define CELLNUM 512000

#define LONG_PUTR(var,dest) \
asm volatile ("putr %0,%1\n"::"r"(var),"r"(dest):"memory")
#define LONG_GETR(var) \
asm volatile ("getr %0\n":"=r"(var)::"memory")
#define LONG_PUTC(var,dest) \
asm volatile ("putc %0,%1\n"::"r"(var),"r"(dest):"memory")
#define LONG_GETC(var) \
asm volatile ("getc %0\n":"=r"(var)::"memory")


__thread_local volatile unsigned long get_reply,put_reply;
//__thread_local volatile unsigned long startEndIdx[2],blockIdx,cellIdx[2];
__thread_local int myId;
//__thread_local int upper_slave[BLOCKNUM64K*BLOCKNUM64K/2],owner_slave[BLOCKNUM64K*BLOCKNUM64K/2],lower_slave[BLOCKNUM64K*BLOCKNUM64K/2],neighbor_slave[BLOCKNUM64K*BLOCKNUM64K/2];

extern SCALAR *upperReordered,*lowerReordered,*diagReordered,*xReordered,*b_MLBSlave;
extern int *ownReordered,*neiReordered,*blockStarts,*cellStarts,*colStarts,*ownNeiSendList;
extern int cpeBlockNum,maxEdges,maxCells,mshBlockNum,spIndex,maxBlock;


void func(){

	int i,j,k,row,col,rowId;
	myId = athread_get_id(-1);
	get_reply = 0;

	volatile int maxXNum = maxBlock;
	volatile int maxCell = maxCells;
	volatile int maxEdge = maxEdges;
	volatile int cpeBlockNums = cpeBlockNum;
	volatile int index = spIndex;
	volatile int totalLength = maxXNum*sizeof(SCALAR)+maxCell*2*sizeof(SCALAR)+maxEdge*4*sizeof(int)+BLOCKNUM64K*6*sizeof(int)+(_total_send_pcg+_total_recv_pcg)*32;
	if(totalLength > 64*1024*0.99 && myId==0){
//	if(myId==0){
		printf("The allocated LDM exceeds 64KB, the memory size is %d\n",totalLength);
	}

	volatile SCALAR b_slave[maxXNum];
	volatile SCALAR diagUpper_slave[maxCell];
	volatile int diagOwner_slave[maxCell];

	volatile int cell_slave[BLOCKNUM64K+1];

	volatile int blockStarts_slave[BLOCKNUM64K*4];
    volatile int startIdx,blockLen,startBlockIdx,endIdx;
	volatile int integral,remainder,kk,l,kkk,jj;
	
	//寄存器通信发送面数据
	volatile SCALAR sLower_slave[maxEdge],rLower_slave[maxEdge];
	volatile int sNeighbor_slave[maxEdge],rNeighbor_slave[maxEdge];
//	if(maxEdge>total_send_pcg*6) printf("%d,%d,%d\n",myId,maxEdge,total_send_pcg);

	volatile int sPackIdx,testNum,nonZeroNum;
	volatile int ownNeiSendIdx[BLOCKNUM64K];


	sPackIdx = 0;
	row = myId;
	row = myId + index*BLOCKNUM64K;
	startBlockIdx = row*(1+2*cpeBlockNums-row)/2;
//	for(i=0;i<maxXNum;i++){b_slave[i]=0;}
	get_reply = 0;
	athread_get(PE_MODE,&blockStarts[4*startBlockIdx],&blockStarts_slave[0],BLOCKNUM64K*4*sizeof(int),&get_reply,0,0,0);
	athread_get(PE_MODE,&cellStarts[index*BLOCKNUM64K],&cell_slave[0],(1+BLOCKNUM64K)*sizeof(int),&get_reply,0,0,0);
//	athread_get(PE_MODE,&ownNeiSendList[row*BLOCKNUM64K],&ownNeiSendList_slave[0],BLOCKNUM64K*sizeof(int),&get_reply,0,0,0);
	while(get_reply!=2);
	get_reply = 0;
	athread_get(PE_MODE,&b_MLBSlave[cell_slave[myId]],&b_slave[0],(cell_slave[myId+1]-cell_slave[myId])*sizeof(SCALAR),&get_reply,0,0,0);
	while(get_reply!=1);


//	if(myId==0) printf("%d,%d\n",cell_slave[1],cell_slave[0]);

	startIdx = blockStarts_slave[6];
	blockLen = myId==BLOCKNUM64K-1? 0 : blockStarts_slave[4*(BLOCKNUM64K-myId-1)+3]-startIdx;
	if(blockLen>0){
		get_reply = 0;
		athread_get(PE_MODE,&upperReordered[startIdx],&sLower_slave[0],blockLen*sizeof(SCALAR),&get_reply,0,0,0);
		athread_get(PE_MODE,&ownReordered[startIdx],&sNeighbor_slave[0],blockLen*sizeof(int),&get_reply,0,0,0);
		while(get_reply!=2);
		for(i=0;i<blockLen;i++){
			b_slave[sNeighbor_slave[i]-cell_slave[myId]] -=sLower_slave[i];
		}
		get_reply = 0;
		athread_get(PE_MODE,&lowerReordered[startIdx],&sLower_slave[0],blockLen*sizeof(SCALAR),&get_reply,0,0,0);	
		athread_get(PE_MODE,&neiReordered[startIdx],&sNeighbor_slave[0],blockLen*sizeof(int),&get_reply,0,0,0);
		while(get_reply!=2);

	}


	get_reply = 0;
	athread_get(PE_MODE,&upperReordered[blockStarts_slave[2]],&diagUpper_slave[0],(blockStarts_slave[3]-blockStarts_slave[2])*sizeof(SCALAR),&get_reply,0,0,0);
	athread_get(PE_MODE,&ownReordered[blockStarts_slave[2]],&diagOwner_slave[0],(blockStarts_slave[3]-blockStarts_slave[2])*sizeof(int),&get_reply,0,0,0);
	while(get_reply!=2);
	for(k=0;k<blockStarts_slave[3]-blockStarts_slave[2];k++){
		b_slave[diagOwner_slave[k]-cell_slave[myId]] -= diagUpper_slave[k];
	}

	get_reply = 0;
	athread_get(PE_MODE,&neiReordered[blockStarts_slave[2]],&diagOwner_slave[0],(blockStarts_slave[3]-blockStarts_slave[2])*sizeof(int),&get_reply,0,0,0);
	athread_get(PE_MODE,&lowerReordered[blockStarts_slave[2]],&diagUpper_slave[0],(blockStarts_slave[3]-blockStarts_slave[2])*sizeof(SCALAR),&get_reply,0,0,0);
	while(get_reply!=2);
	for(k=0;k<blockStarts_slave[3]-blockStarts_slave[2];k++){
		b_slave[diagOwner_slave[k]-cell_slave[myId]] -= diagUpper_slave[k];
	}


	//communicate Lower
	for(j=0;j<BLOCKNUM64K;j++){ownNeiSendIdx[j]=0;}
	for(j=0;j<_total_send_pcg;j++){
		if(myId>_sPacks[j].dst_id) continue;
		startIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+2]+ownNeiSendIdx[_sPacks[j].dst_id];
		endIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+3];
		kk=startIdx-blockStarts_slave[6];
		if(endIdx-startIdx>=6){
			for(i=0;i<6;i++){
				_sPacks[j].data[i]=sLower_slave[kk+i];
			}
			ownNeiSendIdx[_sPacks[j].dst_id]+=6;
			_sPacks[j].indM = 6;
		}else{
			for(i=0;i<endIdx-startIdx;i++){
				_sPacks[j].data[i]=sLower_slave[kk+i];
			}
			ownNeiSendIdx[_sPacks[j].dst_id]+=endIdx-startIdx;
			_sPacks[j].indM = endIdx-startIdx;
		}
	}

	transform_data();
	nonZeroNum=0;
	for(i=0;i<_total_recv_pcg;i++){
		if(_rPacks[i].src_id<myId){
			for(j=0;j<_rPacks[i].indM;j++){
				rLower_slave[nonZeroNum] = _rPacks[i].data[j];
				nonZeroNum++;
			}
		}
	}

	//communicate Neighbor
	for(j=0;j<BLOCKNUM64K;j++){ownNeiSendIdx[j]=0;}
	for(j=0;j<_total_send_pcg;j++){
		if(myId>_sPacks[j].dst_id) continue;
		startIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+2]+ownNeiSendIdx[_sPacks[j].dst_id];
		endIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+3];
		kk=startIdx-blockStarts_slave[6];
		if(endIdx-startIdx>=6){
			for(i=0;i<6;i++){
				_sPacks[j].data[i]=sNeighbor_slave[kk+i];
			}
			ownNeiSendIdx[_sPacks[j].dst_id]+=6;
			_sPacks[j].indM = 6;
		}else{
			for(i=0;i<endIdx-startIdx;i++){
				_sPacks[j].data[i]=sNeighbor_slave[kk+i];
			}
			ownNeiSendIdx[_sPacks[j].dst_id]+=endIdx-startIdx;
			_sPacks[j].indM = endIdx-startIdx;
		}
	}


	transform_data();
	nonZeroNum=0;
	int lastId=0;
	k=0;
	for(i=0;i<_total_recv_pcg;i++){
		if(_rPacks[i].src_id<myId){
			for(j=0;j<_rPacks[i].indM;j++){
				rNeighbor_slave[nonZeroNum] = _rPacks[i].data[j];
				nonZeroNum++;
			}
		}
	}

	//计算稠密块中下三角
	for(i=0;i<nonZeroNum;i++){
		b_slave[rNeighbor_slave[i]-cell_slave[myId]] -= rLower_slave[i];
//		if(rNeighbor_slave[i]==323546) printf("rNeighbor: %f,%f\n",rLower_slave[i],b_slave[rNeighbor_slave[i]-cell_slave[myId]]);
	}

	put_reply = 0;
	athread_put(PE_MODE,&b_slave[0],&b_MLBSlave[cell_slave[myId]],(cell_slave[myId+1]-cell_slave[myId])*sizeof(SCALAR),&put_reply,0,0);
	while(put_reply!=1);
}
