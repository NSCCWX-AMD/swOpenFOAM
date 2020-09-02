#include "swNegSumDiag_host.h"
#include "BlockOrdering.h"
#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

void BlockOrderingMLB(struct LDUMatrix Matrix){

	int i,j,blockIdx,cellNum,edgeNum,rowId,Idx;
	postCellOrder = (int*)malloc(sizeof(int)*Matrix.numCell);
	postEdgeOrder = (int*)malloc(sizeof(int)*Matrix.numEdge);
	printf("cells:%d,edges:%d\n",Matrix.numCell,Matrix.numEdge);
	int totalSize = Matrix.numCell*sizeof(SCALAR)+Matrix.numCell*4*2/3*(sizeof(int)+sizeof(SCALAR))+Matrix.numCell*4/3*2*(sizeof(int)+sizeof(SCALAR))+Matrix.numCell*4/3/(24/sizeof(SCALAR))*32;
	printf("The estimated mshBlockNum = %d\n",totalSize/64/64/1024);
//	if(ceil(((double)(Matrix.numCell+Matrix.numEdge))*4.0/4.0/1024/1024)>2){
//		mshBlockNum = ceil(((double)(Matrix.numCell+Matrix.numEdge))*4.0/4.0/1024/1024);
//	}else{
	mshBlockNum = totalSize/64/64/1024+4;
//	}
	printf("The dataset is divided into %d parts to fullfill the memory restriction of LDM\n",mshBlockNum);
	int* cellWeights = (int*)malloc(sizeof(int)*Matrix.numCell);
	int* edgeWeights = (int*)malloc(sizeof(int)*Matrix.numEdge);
	for(i=0;i<Matrix.numCell;i++) cellWeights[i] = 1;
	for(i=0;i<Matrix.numEdge;i++) edgeWeights[i] = 2;
	MLB_graph graph;
	graph.owner = Matrix.rowAddr;
	graph.neighbor = Matrix.colAddr;
	graph.cellWeights = cellWeights;
	graph.edgeWeights = edgeWeights;
	graph.cellNum = Matrix.numCell;
	graph.edgeNum = Matrix.numEdge;

	cpeBlockNum = mshBlockNum*BLOCKNUM64K;
	mtxBlockNum = cpeBlockNum*(cpeBlockNum+1)/2;
	int levels = 2;
	int *blockNums = (int*)malloc(sizeof(int)*2);
	blockNums[0] = mshBlockNum;
	blockNums[1] = BLOCKNUM64K;
	blockStarts = (int*)malloc(sizeof(int)*mtxBlockNum*4);
	cellStarts = (int*)malloc(sizeof(int)*(cpeBlockNum+1));

	MLB_Multilevel_ordering(graph,levels,blockNums,blockStarts,cellStarts,postCellOrder,postEdgeOrder);

	ownReordered = (int*)malloc(sizeof(int)*Matrix.numEdge);
	neiReordered = (int*)malloc(sizeof(int)*Matrix.numEdge);
	upperReordered = (SCALAR*)malloc(sizeof(SCALAR)*Matrix.numEdge);
	lowerReordered = (SCALAR*)malloc(sizeof(SCALAR)*Matrix.numEdge);
	diagReordered = (SCALAR*)malloc(sizeof(SCALAR)*Matrix.numCell);
	xReordered = (SCALAR*)malloc(sizeof(SCALAR)*Matrix.numCell);
	for(i=0;i<Matrix.numEdge;i++){
		if(postEdgeOrder[i]<0){
			ownReordered[-postEdgeOrder[i]] = postCellOrder[Matrix.colAddr[i]];
			neiReordered[-postEdgeOrder[i]] = postCellOrder[Matrix.rowAddr[i]];
			upperReordered[-postEdgeOrder[i]] = Matrix.lower[i];
			lowerReordered[-postEdgeOrder[i]] = Matrix.upper[i];
		}else{
			ownReordered[postEdgeOrder[i]] = postCellOrder[Matrix.rowAddr[i]];
			neiReordered[postEdgeOrder[i]] = postCellOrder[Matrix.colAddr[i]];
			upperReordered[postEdgeOrder[i]] = Matrix.upper[i];
			lowerReordered[postEdgeOrder[i]] = Matrix.lower[i];
		}
	}
	for(i=0;i<Matrix.numCell;i++){
		diagReordered[postCellOrder[i]] = Matrix.diag[i];
	}

	maxEdges = 0;
	maxCells = 0;
	edgeNum = 0;
	cellNum = 0;
	for(Idx=0;Idx<mshBlockNum;Idx++){
		for(i=0;i<BLOCKNUM64K;i++){
			rowId = Idx*BLOCKNUM64K+i;
			for(j=rowId+1;j<(Idx+1)*BLOCKNUM64K;j++){
				blockIdx = rowId*(1+2*cpeBlockNum-rowId)/2+j-rowId;
				edgeNum += blockStarts[4*blockIdx+3]-blockStarts[4*blockIdx+2];
			}
			blockIdx = rowId*(1+2*cpeBlockNum-rowId)/2;
			cellNum = blockStarts[4*blockIdx+3]-blockStarts[4*blockIdx+2];
			maxCells = maxCells > cellNum ? maxCells : cellNum;
			maxEdges = maxEdges > edgeNum ? maxEdges : edgeNum;
			edgeNum = 0;
		}
	}
	maxBlock = Matrix.numCell/cpeBlockNum*2;
	printf("maxEdges = %d\n",maxEdges);
	printf("maxCells = %d\n",maxCells);
	printf("maxBlock = %d\n",maxBlock);

	free(cellWeights);
	free(edgeWeights);
}


void SMVM_SP_BlockOrder(struct LDUMatrix Matrix){
	//BlockOrdering数据准备
//	BlockOrderingMLB(Matrix,x);
	
	int i,edgeNum,cellNum,j,blockIdx,startIdx,endIdx,k,row,col,idxMP,face;
	struct timeval start,end;
	b_MLBSlave   = (SCALAR*)malloc(sizeof(SCALAR)*Matrix.numCell);
//	b_MLBSlaveMP = (SCALAR*)malloc(sizeof(SCALAR)*Matrix.numCell);
	for(i=0;i<Matrix.numCell;i++){b_MLBSlave[postCellOrder[i]]=Matrix.diag[i];}
	int totalSize = (Matrix.numCell/cpeBlockNum*2*2+maxEdges)*sizeof(SCALAR)+(maxEdges*2+cpeBlockNum+1)*sizeof(int);
	if(totalSize > 64*1024){
		printf("Error: the size of arrays exceeds the limits of LDM!\n TotalSize=%d\n",totalSize);
	}
	for(face=0;face<Matrix.numEdge;face++){
		//printf("%d,%d,%d\n",face,l[face],u[face]);
		//assert(u[face]<=face);
		//assert(l[face]<=face);
		Matrix.diag[Matrix.rowAddr[face]] -= Matrix.lower[face];
		Matrix.diag[Matrix.colAddr[face]] -= Matrix.upper[face];
	}


//	initOwnNeiSendList();


//	athread_init();
//	SMVM_reg(ownNeiSendList,mshBlockNum);
//	athread_spawn(initTable,schedule_data);
	gettimeofday(&start,NULL);
	for(spIndex=0;spIndex<mshBlockNum;spIndex++){
//		printf("init:%d\n",spIndex);
		initTable(spIndex);
//		printf("transform:%d\n",spIndex);
		athread_spawn(func,0);
		idxMP = spIndex==mshBlockNum-1? 0:spIndex+1;
		for(j=0;j<BLOCKNUM64K;j++){
			row = idxMP*BLOCKNUM64K+j;
			col=(idxMP+1)*BLOCKNUM64K;
			if(col>=cpeBlockNum) break;
			blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
			startIdx = blockStarts[4*blockIdx+2];
			col=cpeBlockNum-1;
			blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
			endIdx   = blockStarts[4*blockIdx+3];
			for(k=startIdx;k<endIdx;k++){
				b_MLBSlave[neiReordered[k]] -= lowerReordered[k];
				b_MLBSlave[ownReordered[k]] -= upperReordered[k];
			}
		}

		athread_join();
//		printf("destroy:%d\n",spIndex);
		destroyTable(spIndex);
	}
//	athread_join();
//	athread_halt();
	gettimeofday(&end,NULL);
	int timeuse = 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
	printf("Slave Processor Time: %f us\n",(double)timeuse);

//	static int flag=0;
//	if(flag==0){
//	for(i=0;i<Matrix.numCell;i++){
//		if(Matrix.diag[i] != b_MLBSlave  [postCellOrder[i]]){
//			printf("%d,%d,%f,%f\n",i,postCellOrder[i],Matrix.diag[i],b_MLBSlave[postCellOrder[i]]);
//		}
//	}
//	flag++;
//	}
//	free(postCellOrder);
//	free(postEdgeOrder);
//	free(blockStarts);
//	free(cellStarts);
//	free(ownReordered);
//	free(neiReordered);
//	free(upperReordered);
//	free(lowerReordered);
//	free(diagReordered);
//	free(xReordered);
	free(b_MLBSlave);

}

void initOwnNeiSendList(){
	int num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	int startIdx,endIdx,startListIdx,i,j,k,row,col,blockIdx,upperIdx,lowerIdx;
	ownNeiSendList = (int*)malloc(sizeof(int)*num);
	for(k=0;k<num;k++){ownNeiSendList[k]=0;}
	for(i=0;i<mshBlockNum;i++){
		for(j=0;j<BLOCKNUM64K;j++){
			row = i*BLOCKNUM64K+j;
			col = row+1;
			while(col < (i+1)*BLOCKNUM64K && col > i*BLOCKNUM64K){
				blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
				startIdx = blockStarts[4*blockIdx+2];
				endIdx   = blockStarts[4*blockIdx+3];
				upperIdx = row * BLOCKNUM64K + col - i * BLOCKNUM64K;
				lowerIdx = col * BLOCKNUM64K + row - i * BLOCKNUM64K;
//				if(row==1 && col==2){
//					printf("%d,%d\n",startIdx,endIdx);
//				}
				col++;
				if(endIdx-startIdx==0) continue;
				ownNeiSendList[upperIdx] = (endIdx - startIdx-1)/6+1;
				ownNeiSendList[lowerIdx] = (endIdx - startIdx-1)/6+1;
				
//				if(upperIdx == 26 || lowerIdx ==1664){
//				  printf("%d\n",row);
//				  printf("%d\n",col);
//				}
//				for(k=startIdx;k<endIdx;k++){
//					b_MLB[ownReordered[k]] += upperReordered[k]*xReordered[neiReordered[k]];
//					if(ownReordered[k]==MONITORIDX)
//					  printf("upper: %f,%f,%f\n",upperReordered[k],xReordered[neiReordered[k]],b_MLB[MONITORIDX]);
//				}
			}
		}
	}
	//check
	for(i=0;i<mshBlockNum;i++){
		for(j=0;j<BLOCKNUM64K;j++){
			row = i*BLOCKNUM64K+j;
			for(col=0;col<BLOCKNUM64K;col++){
				upperIdx = row * BLOCKNUM64K + col;
				lowerIdx = (col + (i * BLOCKNUM64K)) * BLOCKNUM64K + row - (i*BLOCKNUM64K);
				if(ownNeiSendList[upperIdx] != ownNeiSendList[lowerIdx]){
					printf("%d,%d,%d,%d\n",row,col,ownNeiSendList[upperIdx],ownNeiSendList[lowerIdx]);
				}
//				if(row==1){
//					printf("%d,%d,%d,%d\n",row,col,ownNeiSendList[upperIdx],ownNeiSendList[lowerIdx]);
//				}

			}
		}
	}


}
