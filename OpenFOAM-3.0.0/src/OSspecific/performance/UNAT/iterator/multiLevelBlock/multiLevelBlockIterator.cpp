#include "multiLevelBlockIterator.hpp"
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <sys/time.h>
#include <fstream>
#include "multiLevelBlockIterator.h"

#ifdef __cplusplus
extern "C"{
#endif
#include "BlockOrdering.h"
#include <athread.h>
void slave_multiLevelBlockIterator_e2v_slave(MLB_edge2VertexPara *para);
void slave_multiLevelBlockIterator_v2e_slave(MLB_vertex2EdgePara *para);
void slave_multiLevelBlockIterator_array_slave(MLB_arrayPara *para);
//#include "edge2VertexIter_host.h"
//#include "vertex2EdgeIter_host.h"
#ifdef __cplusplus
}
#endif

//#define CELLDIMS 2
//#define EDGEDIMS 3
#define MAXLDM (50000)
#define MINLDM (46000)
#define MAXMLBITER 5

namespace UNAT
{

// ***********************Constructor*************************************
MultiLevelBlockIterator::MultiLevelBlockIterator(Topology &topo,
			swInt* vertexWeights, swInt* edgeWeights, bool duplicate)
		: Iterator(topo, vertexWeights, edgeWeights, duplicate)
{
//	LOG("MultiLevelBlockIterator");
    //printf("start MultiLevelBlocking reorder...\n");
	MLBReorder(topo, 2, 0, 0);
	//printf("The domain is divided into %d*%d blocks.\n",
	//			this->_mshBlockNum,this->_mshBlockNum);

	swInt vertexNumber = this->getTopology()->getVertexNumber();
	swInt edgeNumber   = this->getTopology()->getEdgeNumber();

	for(int i=0;i<edgeNumber;i++)
	{
		this->getEdgeMap()
			.insert(pair<swInt,swInt>(i,this->_postEdgeOrder[i]));
	}
	for(int i=0;i<vertexNumber;i++)
	{
		this->getVertexMap()
			.insert(pair<swInt,swInt>(i,this->_postVertexOrder[i]));
	}
	printf("inner topology\n");
	this->reformInnerTopology();
	this->initOwnNeiSendList();

//	printArray("%d",this->_blockStarts,this->_mtxBlockNum);
}

void MultiLevelBlockIterator::MLBReorder(Topology &topo, swInt ref,
			swInt lastBlockNum, swInt iter)
{
	swInt* vertexWeights;
	swInt* edgeWeights;
	MLB_graph graph;
	swInt* blockNums;
	swInt  levels;
	swInt vertexNumber = this->getTopology()->getVertexNumber();
	swInt edgeNumber   = this->getTopology()->getEdgeNumber();

	this->_postVertexOrder = (swInt*)malloc(sizeof(swInt)*vertexNumber);	
	this->_postEdgeOrder = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	swInt* postVertexOrder = this->_postVertexOrder;;
	swInt* postEdgeOrder = this->_postEdgeOrder;;
	swInt totalSize
		// sendX, recvX, b
		= vertexNumber*2*sizeof(swFloat)
		// 对角块的owner，neighbour和非零元素
		+ edgeNumber*2*2*(sizeof(swInt)+sizeof(swFloat))
		// 发送接收计数数组和blockStart，cellStarts
		+ BLOCKNUM64K*BLOCKNUM64K*5*sizeof(swInt);

	this->_mshBlockNum = MAX(totalSize/64/64/1024+1,ref);
	this->_mshBlockNum = ref;
	//printf("The estimated mshBlockNum = %d\n",this->_mshBlockNum);

	vertexWeights = (swInt*)malloc(sizeof(swInt)*vertexNumber);
	edgeWeights = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	for(int i=0;i<edgeNumber;i++) edgeWeights[i] = 2;
	for(int i=0;i<vertexNumber;i++) vertexWeights[i] = 1;

	graph.owner       = this->getTopology()->getStartVertices();
	graph.neighbor    = this->getTopology()->getEndVertices();
	graph.cellWeights = this->getVertexWeights();
	graph.edgeWeights = this->getEdgeWeights();
	graph.cellNum     = vertexNumber;
	graph.edgeNum     = edgeNumber;

	this->_cpeBlockNum = this->_mshBlockNum*BLOCKNUM64K;
	this->_mtxBlockNum = this->_cpeBlockNum*(this->_cpeBlockNum+1)/2;
	levels = 2;
	blockNums = (swInt*)malloc(sizeof(swInt)*levels);
	blockNums[0] = this->_mshBlockNum;
	blockNums[1] = BLOCKNUM64K;
	this->_blockStarts = (swInt*)malloc(sizeof(swInt)*this->_mtxBlockNum*4);
	this->_vertexStarts =
		(swInt*)malloc(sizeof(swInt)*(this->_cpeBlockNum+1));

	MLB_Multilevel_ordering(graph,levels,blockNums,this->_blockStarts,
				this->_vertexStarts,postVertexOrder,postEdgeOrder);

	// Get the receive counts of each segment
	swInt *blockStarts = this->_blockStarts;
	swInt maxRecvNum = 0;
	this->_recvStarts = (swInt*)malloc(this->_cpeBlockNum*sizeof(swInt));
	swInt *recvStarts = this->_recvStarts;
	for(int iseg=0;iseg<this->_cpeBlockNum;iseg++)
	{
		recvStarts[iseg] = 0;
	}
	for(int irow=0;irow<this->_cpeBlockNum;irow++)
	{
		int iseg = irow/BLOCKNUM64K;
		int localRow = irow -iseg*BLOCKNUM64K;
		for(int icol=localRow+1;icol<BLOCKNUM64K;icol++)
		{
			int blockIdx
				= irow*(1+2*this->_cpeBlockNum-irow)/2
				+ icol+iseg*BLOCKNUM64K-irow;
			recvStarts[icol+iseg*BLOCKNUM64K]
				+= blockStarts[4*blockIdx+3]-blockStarts[4*blockIdx+2];
		}
	}

	this->_maxEdges = 0;
	this->_maxCells = 0;
	this->_maxXNum  = 0;
	swInt xNum    = 0;
	swInt edgeNum = 0;
	swInt cellNum = 0;
	swInt recvNum = 0;
	swInt blockIdx;
	for(int i=0;i<this->_cpeBlockNum;i++)
	{
		for(int j=i+1;j<this->_cpeBlockNum;j++)
		{
			blockIdx = i*(1+2*this->_cpeBlockNum-i)/2+j-i;
			edgeNum += this->_blockStarts[4*blockIdx+3]
				- this->_blockStarts[4*blockIdx+2];
		}
		blockIdx = i*(1+2*this->_cpeBlockNum-i)/2;
		cellNum = this->_blockStarts[4*blockIdx+3]
			- this->_blockStarts[4*blockIdx+2];
		xNum = this->_vertexStarts[i+1]-this->_vertexStarts[i];
		recvNum = this->_recvStarts[i];
		int edgeNums = edgeNum + cellNum;
		int cellNums = edgeNum + xNum;
		this->_maxXNum
			= this->_maxXNum > xNum ? this->_maxXNum : xNum;	
		this->_maxCells
			= this->_maxCells > cellNums ? this->_maxCells : cellNums;
		this->_maxEdges
			= this->_maxEdges > edgeNums ? this->_maxEdges : edgeNums;
		edgeNum = 0;
	}
	this->_maxEdgesUnsymm = (this->_maxCells*2+this->_maxEdges)*1.01;
//printf("maxXNum: %d,maxCells: %d,maxEdges: %d\n",this->_maxXNum,this->_maxCells,this->_maxEdges);
	int cellDims = this->getVertexWeights()[0];
	int edgeDims = this->getEdgeWeights()[0];
	swInt totalLength
		= this->_maxCells*cellDims*sizeof(swFloat)
//		+ this->_maxCells*sizeof(siwFloat)
		+ this->_maxEdges*2*sizeof(swInt)
		+ this->_maxEdges*edgeDims*sizeof(swFloat)
		+ BLOCKNUM64K*5*sizeof(swInt);
	printf("%dth TotalLength: %d, blockNum: %d\n", iter, totalLength, this->_mshBlockNum);
	iter++;
	if(totalLength > MAXLDM)
	{
		swInt newBlockNum = (totalLength*this->_mshBlockNum-1)/MAXLDM+1;
		if( newBlockNum <= this->_mshBlockNum) newBlockNum = this->_mshBlockNum + 1;
		//printf("\t%d parts increased to %d parts\n", this->_mshBlockNum, newBlockNum);
	   	MLBReorder(topo, newBlockNum, this->_mshBlockNum, iter);
	}else if(totalLength < MINLDM && this->_mshBlockNum>2 && iter <= MAXMLBITER)
	{
		swInt newBlockNum = (totalLength*this->_mshBlockNum)/MINLDM;
		if( newBlockNum >= this->_mshBlockNum) newBlockNum = this->_mshBlockNum - 1;
		//printf("\t%d parts decreased to %d parts\n", this->_mshBlockNum, newBlockNum);
	   	MLBReorder(topo, MAX(2, newBlockNum ), this->_mshBlockNum, iter);
	}

}

void MultiLevelBlockIterator::reorderEdges(swInt* startVertices,
			swInt* endVertices, swInt edgeNumber, swInt vertexNumber)
{
	if(edgeNumber!=this->getTopology()->getEdgeNumber())
	{
		LOG("The edge arrays do not match the topology!");
	}

	map<swInt, swInt>::iterator iter;
	this->_owner   = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	this->_neighbor = (swInt*)malloc(sizeof(swInt)*edgeNumber);

	for(iter = this->getEdgeMap().begin();
				iter!=this->getEdgeMap().end();iter++)
	{
//if(iter->first==0) printf("%d\n",iter->second);
		if(iter->second<0)
		{
			this->_owner[-iter->second-1]
				= this->getVertexMap()[endVertices[iter->first]];
			this->_neighbor[-iter->second-1]
				= this->getVertexMap()[startVertices[iter->first]];
		} else
		{
			this->_owner[iter->second-1]
				= this->getVertexMap()[startVertices[iter->first]];
			this->_neighbor[iter->second-1]
				= this->getVertexMap()[endVertices[iter->first]];
		}
	}

	for(int i=0;i<edgeNumber;i++)
	{
		startVertices[i] = this->_owner[i];
		endVertices[i]   = this->_neighbor[i];
	}

//	this->writeTopology();
	
//	LOG("reorderEdgesFromEdge");
}

void MultiLevelBlockIterator::reorderNeighbor(
			swInt* firstEdgeVertices, swInt* vertexNeighbours,
			swInt edgeNumber, swInt vertexNumber)
{
	this->_firstEdgeVertices
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	this->_vertexNeighbours
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	for(int i=0;i<edgeNumber*2;i++) 
	{
		this->_firstEdgeVertices[i] = firstEdgeVertices[i];
		this->_vertexNeighbours[i]  = vertexNeighbours[i];
	}

	swInt cpeBlockNum   = this->getCpeBlockNum();
	swInt mshBlockNum   = this->getMshBlockNum();
	swInt mtxBlockNum   = this->getMtxBlockNum();
	swInt* blockStarts  = this->getBlockStarts();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getTopology()->getStartVertices();
	swInt* neighbor     = this->getTopology()->getEndVertices();

	swInt* vertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* accuVertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* tmpVN
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	swInt* tmpFEV
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	this->_blockStartsUnsymm
		= (swInt*)malloc(sizeof(swInt)*cpeBlockNum*cpeBlockNum*4);
	swInt* blockStartsUnsymm = this->_blockStartsUnsymm;

	int blockIdx, startBlockIdx, endBlockIdx, row1, col1, blockIdxUnsymm;
	int accuEdges= 0;
	int startIdx = 0;
	int endIdx   = 0;
	for(int i=0;i<mshBlockNum;i++)
	{
		for(int j=0;j<BLOCKNUM64K;j++)
		{
			int row = i*BLOCKNUM64K+j;
			int cellLen = vertexStarts[row+1]-vertexStarts[row];
			assert(cellLen<vertexNumber/cpeBlockNum*2);
			for(int col=0;col<cpeBlockNum;col++)
			{
				blockIdxUnsymm = row*cpeBlockNum+col;
				blockStartsUnsymm[4*blockIdxUnsymm  ] = row;	
				blockStartsUnsymm[4*blockIdxUnsymm+1] = col;	
				for(int k=0;k<vertexNumber/cpeBlockNum*2;k++)
				{
					vertexEdgeNumbers[k]=0;
					accuVertexEdgeNumbers[k]=0;
				}
				if(col==row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
		    			vertexEdgeNumbers[owner[k]-vertexStarts[row]]++;
		    			vertexEdgeNumbers[neighbor[k]-vertexStarts[row]]++;
		    		}
					for(int k=1;k<cellLen;k++)
					{
						assert(vertexEdgeNumbers[k-1]>=0);
						accuVertexEdgeNumbers[k+1]
							= accuVertexEdgeNumbers[k]
							+ vertexEdgeNumbers[k-1];
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = neighbor[k]+1-vertexStarts[row];
						col1 = owner[k];
						tmpVN[accuVertexEdgeNumbers[row1]+startIdx]=col1;
						tmpFEV[accuVertexEdgeNumbers[row1]+startIdx]
							=neighbor[k];
						accuVertexEdgeNumbers[row1]++;
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = owner[k]+1-vertexStarts[row];
						col1 = neighbor[k];
						tmpVN[accuVertexEdgeNumbers[row1]+startIdx]=col1;
						tmpFEV[accuVertexEdgeNumbers[row1]+startIdx]
							=owner[k];
						accuVertexEdgeNumbers[row1]++;
					}
					blockStartsUnsymm[4*blockIdxUnsymm+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx)*2;
					blockStartsUnsymm[4*blockIdxUnsymm+3] = startIdx;	
				} else if(col < row)
				{
     				blockIdx = col*(1+2*this->_cpeBlockNum-col)/2+row-col;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
						tmpVN[k-startBlockIdx+startIdx] = owner[k];
						tmpFEV[k-startBlockIdx+startIdx] = neighbor[k];
		    		}
					blockStartsUnsymm[4*blockIdxUnsymm+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx);
					blockStartsUnsymm[4*blockIdxUnsymm+3] = startIdx;	
				} else if(col > row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2+col-row;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
						tmpVN[k-startBlockIdx+startIdx] = neighbor[k];
						tmpFEV[k-startBlockIdx+startIdx] = owner[k];
		    		}
					blockStartsUnsymm[4*blockIdxUnsymm+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx);
					blockStartsUnsymm[4*blockIdxUnsymm+3] = startIdx;	
				}
			}
		}
	}
	for(int i=0;i<edgeNumber*2;i++)
	{
		firstEdgeVertices[i] = tmpFEV[i];
		vertexNeighbours[i] = tmpVN[i];
	}
	swInt maxEdgesUnsymm = 0;
	swInt edgeNum = 0;
	for(int i=0;i<cpeBlockNum;i++)
	{
		edgeNum = blockStartsUnsymm[4*((i+1)*cpeBlockNum-1)+3]
			- blockStartsUnsymm[4*i*cpeBlockNum+2];
		maxEdgesUnsymm = maxEdgesUnsymm>edgeNum?maxEdgesUnsymm:edgeNum;
	}
	assert(this->_maxEdgesUnsymm>maxEdgesUnsymm);
	free(tmpFEV);
	free(tmpVN);
	free(vertexEdgeNumbers);
	free(accuVertexEdgeNumbers);
//	assert(startIdx==edgeNumber*2);
//	for(int i=0;i<edgeNumber*2;i++)
//	{
//		if(vertexNeighbours[i]==-1 || firstEdgeVertices[i]==-1)
//		{
//			printf("%d,%d,%d\n",i,vertexNeighbours[i],firstEdgeVertices[i]);
//			exit(1);
//		}
//	}
//
//	for(int i=0;i<mshBlockNum;i++)
//	{
//		for(int j=0;j<BLOCKNUM64K;j++)
//		{
//			int row = i*BLOCKNUM64K+j;
//			for(int col=0;col<cpeBlockNum;col++)
//			{
//				int blockIdxV = row*cpeBlockNum+col;
//				printf("(%4d,%4d),",blockStartsUnsymm[4*blockIdxV+2],blockStartsUnsymm[4*blockIdxV+3]);
//			}
//		}
//	}
//printf("row:%d, col:%d, edgeNums:%d, vertexNums:%d, startIdx:%d\n",row,col,endBlockIdx-startBlockIdx,cellLen,startIdx);
//for(int m=0;m<140;m++)
//{
//	for(int n=0;n<15;n++)
//	{
//		printf("%4d ",owner[m*15+n]);
//	}
//	printf("\n");
//}
//printf("********************************************\n");
	LOG("reorderEdgesFromVertex");
}

void MultiLevelBlockIterator::arrayIteration(Arrays* paraData,
			coupledArrayOperator* cOpt, int optNum)
{
	MLB_arrayPara para =
	{
		this->_vertexStarts,
		this->_blockStarts,
		this->getTopology()->getVertexNumber(),
		this->getTopology()->getEdgeNumber(),
		this->_mshBlockNum,

		1,

		paraData,

		cOpt,
		optNum
	};
//	for(int spIndex=this->_mshBlockNum-1;spIndex>=0;spIndex--)
//	{
//		para.spIndex = spIndex;
		__real_athread_spawn((void*)slave_multiLevelBlockIterator_array_slave, &para);
		athread_join();
//	}
}

void MultiLevelBlockIterator::edge2VertexIteration(Arrays* paraData,
			coupledOperator* cOpt, int optNum)
{
//	double time1,time2;
//	fun_host(backEdgeData, frontEdgeData, selfConnData, vertexData,
//				this->getTopology()->getStartVertices(),
//				this->getTopology()->getEndVertices());
	MLB_edge2VertexPara para = 
	{
		this->getTopology()->getStartVertices(),
		this->getTopology()->getEndVertices(),
		this->_vertexStarts,
		this->_blockStarts,
		this->_recvStarts,
		this->getTopology()->getVertexNumber(),
		this->getTopology()->getEdgeNumber(),
		this->_mshBlockNum,
		this->_maxEdges+this->_maxCells,
		this->_maxXNum,

		1,

		this->_schedule_data,

		paraData,

		cOpt,
		optNum
	};

//	getTime(time1);
	swInt cpeBlockNum  = this->_cpeBlockNum;
	swInt mshBlockNum  = this->_mshBlockNum;
	swInt *blockStarts = this->_blockStarts;
	swInt *owner       = this->getTopology()->getStartVertices();
	swInt *neighbor    = this->getTopology()->getEndVertices();

	Arrays backEdgeData_m[optNum], frontEdgeData_m[optNum];
	Arrays selfConnData_m[optNum], vertexData_m[optNum];
	Arrays backEdgeData_p[optNum],frontEdgeData_p[optNum];
	Arrays selfConnData_p[optNum], vertexData_p[optNum];
	FieldData data[optNum];
	Arrays *backEdgeArrays[optNum], *frontEdgeArrays[optNum];
	Arrays *backEdgeArrays_p[optNum], *frontEdgeArrays_p[optNum];
	for(int iOpt=0;iOpt<optNum;iOpt++)
	{
		constructSingleArrays(&backEdgeData_m[iOpt],
					cOpt[iOpt].data->backEdgeData);
		constructSingleArrays(&frontEdgeData_m[iOpt],
					cOpt[iOpt].data->frontEdgeData);
		constructSingleArrays(&selfConnData_m[iOpt],
					cOpt[iOpt].data->selfConnData);
		constructSingleArrays(&vertexData_m[iOpt],
					cOpt[iOpt].data->vertexData);
		data[iOpt].backEdgeData  = &backEdgeData_m[iOpt];
		data[iOpt].frontEdgeData = &frontEdgeData_m[iOpt];
		data[iOpt].selfConnData  = &selfConnData_m[iOpt];
		data[iOpt].vertexData    = &vertexData_m[iOpt];
		backEdgeArrays[iOpt]  = cOpt[iOpt].data->backEdgeData;
		frontEdgeArrays[iOpt] = cOpt[iOpt].data->frontEdgeData;

		constructSingleArrays(&backEdgeData_p[iOpt],
					cOpt[iOpt].data_p->backEdgeData);
		constructSingleArrays(&frontEdgeData_p[iOpt],
					cOpt[iOpt].data_p->frontEdgeData);
		constructSingleArrays(&selfConnData_p[iOpt],
					cOpt[iOpt].data_p->selfConnData);
		constructSingleArrays(&vertexData_p[iOpt],
					cOpt[iOpt].data_p->vertexData);
		backEdgeArrays_p[iOpt]  = cOpt[iOpt].data_p->backEdgeData;
		frontEdgeArrays_p[iOpt] = cOpt[iOpt].data_p->frontEdgeData;
	}

	for(int spIndex=this->_mshBlockNum-1;spIndex>=0;spIndex--)
	{
		para.spIndex = spIndex;
		__real_athread_spawn((void*)slave_multiLevelBlockIterator_e2v_slave, &para);

		int mpIndex = spIndex+1;
		if(mpIndex<mshBlockNum)
		{
			swInt col = (mpIndex+1)*BLOCKNUM64K;
			if(col<cpeBlockNum) 
			{
				int j;
				for(j=0;j<BLOCKNUM64K;j++)
				{
					col = (mpIndex+1)*BLOCKNUM64K;
					swInt row = mpIndex*BLOCKNUM64K+j;
					swInt blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
					swInt startIdx = blockStarts[4*blockIdx+2];
					col = cpeBlockNum-1;
					blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
					swInt endIdx = blockStarts[4*blockIdx+3];

					for(int iOpt=0;iOpt<optNum;iOpt++)
					{
						backEdgeData_m[iOpt].fArraySizes = endIdx-startIdx;
						frontEdgeData_m[iOpt].fArraySizes= endIdx-startIdx;
						selfConnData_m[iOpt].fArraySizes = 0;
						backEdgeData_p[iOpt].fArraySizes = endIdx-startIdx;
						frontEdgeData_p[iOpt].fArraySizes =endIdx-startIdx;

						for(int iArray=0;iArray<backEdgeData_m[iOpt].fArrayNum;iArray++)
						{
							swInt dims=getArrayDims(&backEdgeData_m[iOpt],
										iArray);
							backEdgeData_m[iOpt].floatArrays[iArray]
								= backEdgeArrays[iOpt]->floatArrays[iArray]
								+ startIdx*dims;
						}
						for(int iArray=0;iArray<frontEdgeData_m[iOpt].fArrayNum;iArray++)
						{
							swInt dims=getArrayDims(&frontEdgeData_m[iOpt],
										iArray);
							frontEdgeData_m[iOpt].floatArrays[iArray]
								=frontEdgeArrays[iOpt]->floatArrays[iArray]
								+startIdx*dims;
						}
						for(int iArray=0;iArray<backEdgeData_p[iOpt].fArrayNum;
									iArray++)
						{
							swInt dims =getArrayDims(&backEdgeData_p[iOpt],
										iArray);
							backEdgeData_p[iOpt].floatArrays[iArray]
							= backEdgeArrays_p[iOpt]->floatArrays[iArray]
								+ startIdx*dims;
						}
						for(int iArray=0;iArray<frontEdgeData_p[iOpt].fArrayNum;
									iArray++)
						{
							swInt dims=getArrayDims(&frontEdgeData_p[iOpt],
										iArray);
							frontEdgeData_p[iOpt].floatArrays[iArray]
							= frontEdgeArrays_p[iOpt]->floatArrays[iArray]
								+ startIdx*dims;
						}
						cOpt[iOpt].fun_host(&backEdgeData_p[iOpt],
									&frontEdgeData_p[iOpt],
									&selfConnData_p[iOpt], 
									&vertexData_p[iOpt], paraData,
									&owner[startIdx], &neighbor[startIdx],
									&data[iOpt]);
					}
				}
			}
		}
		athread_join();
	}

	int j;
	for(j=0;j<BLOCKNUM64K;j++)
	{
		int col = BLOCKNUM64K;
		swInt row = j;
		swInt blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
		swInt startIdx = blockStarts[4*blockIdx+2];
		col = cpeBlockNum-1;
		blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
		swInt endIdx = blockStarts[4*blockIdx+3];
		for(int iOpt=0;iOpt<optNum;iOpt++)
		{
			backEdgeData_m[iOpt].fArraySizes = endIdx-startIdx;
			frontEdgeData_m[iOpt].fArraySizes = endIdx-startIdx;
			selfConnData_m[iOpt].fArraySizes = 0;
			backEdgeData_p[iOpt].fArraySizes = endIdx-startIdx;
			frontEdgeData_p[iOpt].fArraySizes = endIdx-startIdx;

			for(int iArray=0;iArray<backEdgeData_m[iOpt].fArrayNum;
						iArray++)
			{
				swInt dims = getArrayDims(&backEdgeData_m[iOpt],iArray);
				backEdgeData_m[iOpt].floatArrays[iArray]
					= backEdgeArrays[iOpt]->floatArrays[iArray]
					+ startIdx*dims;
			}
			for(int iArray=0;iArray<frontEdgeData_m[iOpt].fArrayNum;
						iArray++)
			{
				swInt dims = getArrayDims(&frontEdgeData_m[iOpt],iArray);
				frontEdgeData_m[iOpt].floatArrays[iArray]
					= frontEdgeArrays[iOpt]->floatArrays[iArray]
					+ startIdx*dims;
			}
			for(int iArray=0;iArray<backEdgeData_p[iOpt].fArrayNum;
						iArray++)
			{
				swInt dims = getArrayDims(&backEdgeData_p[iOpt],iArray);
				backEdgeData_p[iOpt].floatArrays[iArray]
					= backEdgeArrays_p[iOpt]->floatArrays[iArray]
					+ startIdx*dims;
			}
			for(int iArray=0;iArray<frontEdgeData_p[iOpt].fArrayNum;
						iArray++)
			{
				swInt dims = getArrayDims(&frontEdgeData_p[iOpt],iArray);
				frontEdgeData_p[iOpt].floatArrays[iArray]
					= frontEdgeArrays_p[iOpt]->floatArrays[iArray]
					+ startIdx*dims;
			}
			cOpt[iOpt].fun_host(&backEdgeData_p[iOpt], &frontEdgeData_p[iOpt],
						&selfConnData_p[iOpt], &vertexData_p[iOpt], 
						paraData, &owner[startIdx],&neighbor[startIdx],
						&data[iOpt]);
		}
	}

//	restoreVertexData(vertexData);
//	restoreEdgeData(backEdgeData,frontEdgeData);
}

void MultiLevelBlockIterator::restoreVertexData(Arrays* vertexData)
{
	map<swInt, swInt>::iterator iter;
	swFloat *tmp;
	swInt vertexNum = vertexData->fArraySizes;
	for(int i=0;i<vertexData->fArrayNum;i++)
	{
		swInt dims = vertexData->fArrayDims[i];
		tmp = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat));
		for(iter = this->getVertexMap().begin();
					iter!=this->getVertexMap().end();iter++)
		{
//if(iter->first==90033) printf("%d\n",iter->second);
			for(int j=0;j<dims;j++)
			{
				tmp[iter->first*dims+j]
					= vertexData->floatArrays[i][iter->second*dims+j];
			}
		}
		for(int ivertex=0;ivertex<vertexNum;ivertex++)
		{
			for(int j=0;j<dims;j++)
			{
				vertexData->floatArrays[i][ivertex*dims+j]
					= tmp[ivertex*dims+j];
			}
		}
		free(tmp);
	}
}

void MultiLevelBlockIterator::vertex2EdgeIteration(Arrays *paraData,
			coupledOperator *cOpt, swInt optNum)
{
#define moveArrayPtr(dstArray, srcArray, steps) \
{ \
	for(int iArray=0;iArray<dstArray[iOpt].fArrayNum;iArray++) \
	{ \
		swInt dims=getArrayDims(&dstArray[iOpt], iArray); \
		dstArray[iOpt].floatArrays[iArray]\
			= srcArray[iOpt]->floatArrays[iArray] + (steps)*dims; \
	} \
}

	MLB_vertex2EdgePara para = 
	{
		this->getTopology()->getFirstEdgeVertices(),
		this->getTopology()->getVertexNeighbours(),
		this->_vertexStarts,
		this->_blockStartsUnsymm,
		this->getTopology()->getEdgeNumber(),
		this->getTopology()->getVertexNumber(),
		this->_mshBlockNum,
		this->_maxEdges,
		this->_maxCells,

		0,

		this->_schedule_data,

		paraData,

		cOpt,
		optNum
	};

	swInt cpeBlockNum  = this->_cpeBlockNum;
	swInt mshBlockNum  = this->_mshBlockNum;
	swInt *blockStarts = this->_blockStartsUnsymm;
	swInt *firstEdgeVertice = this->getTopology()->getFirstEdgeVertices();
	swInt *vertexNeighbor   = this->getTopology()->getVertexNeighbours();

	Arrays backEdgeData_m[optNum], frontEdgeData_m[optNum];
	Arrays selfConnData_m[optNum], vertexData_m[optNum];
	FieldData data[optNum];
	Arrays *backEdgeArrays[optNum], *frontEdgeArrays[optNum];
	for(int iOpt=0;iOpt<optNum;iOpt++)
	{
		constructSingleArrays(&backEdgeData_m[iOpt],
					cOpt[iOpt].data->backEdgeData);
		constructSingleArrays(&frontEdgeData_m[iOpt],
					cOpt[iOpt].data->frontEdgeData);
		constructSingleArrays(&selfConnData_m[iOpt],
					cOpt[iOpt].data->selfConnData);
		constructSingleArrays(&vertexData_m[iOpt],
					cOpt[iOpt].data->vertexData);
		data[iOpt].backEdgeData  = &backEdgeData_m[iOpt];
		data[iOpt].frontEdgeData = &frontEdgeData_m[iOpt];
		data[iOpt].selfConnData  = &selfConnData_m[iOpt];
		data[iOpt].vertexData    = &vertexData_m[iOpt];
		backEdgeArrays[iOpt]  = cOpt[iOpt].data->backEdgeData;
		frontEdgeArrays[iOpt] = cOpt[iOpt].data->frontEdgeData;
	}

	swInt row, col, blockIdx, startIdx, endIdx, i;
	for(int spIndex=this->_mshBlockNum-1;spIndex>=0;spIndex--)
	{
		para.spIndex = spIndex;
		__real_athread_spawn((void*)slave_multiLevelBlockIterator_v2e_slave, 
					&para);
		int mpIndex = spIndex+1;
		if(mpIndex<mshBlockNum)
		{
			for(i=0;i<BLOCKNUM64K;i++)
			{
				// lower triangle
				row = mpIndex*BLOCKNUM64K+i;
				col = 0;
				blockIdx = row*cpeBlockNum+col;
				startIdx = blockStarts[4*blockIdx+2];
				col = mpIndex*BLOCKNUM64K;
				blockIdx = row*cpeBlockNum+col;
				endIdx = blockStarts[4*blockIdx+2];

				for(int iOpt=0;iOpt<optNum;iOpt++)
				{
					backEdgeData_m[iOpt].fArraySizes = 0;
					frontEdgeData_m[iOpt].fArraySizes= endIdx-startIdx;
					selfConnData_m[iOpt].fArraySizes = 0;

					moveArrayPtr(frontEdgeData_m, frontEdgeArrays, startIdx);
					cOpt[iOpt].fun_host(&backEdgeData_m[iOpt], &frontEdgeData_m[iOpt], &selfConnData_m[iOpt], &vertexData_m[iOpt], paraData,
								&firstEdgeVertice[startIdx], &vertexNeighbor[startIdx],
								&data[iOpt]);
				}
				// upper triangle
				col = (mpIndex+1)*BLOCKNUM64K;
				if(col>=cpeBlockNum) continue;
				blockIdx = row*cpeBlockNum+col;
				startIdx = blockStarts[4*blockIdx+2];
				col = cpeBlockNum-1;
				blockIdx = row*cpeBlockNum+col;
				endIdx = blockStarts[4*blockIdx+3];

				for(int iOpt=0;iOpt<optNum;iOpt++)
				{
					backEdgeData_m[iOpt].fArraySizes = 0;
					frontEdgeData_m[iOpt].fArraySizes= endIdx-startIdx;
					selfConnData_m[iOpt].fArraySizes = 0;

					moveArrayPtr(frontEdgeData_m, frontEdgeArrays, startIdx);
					cOpt[iOpt].fun_host(&backEdgeData_m[iOpt], &frontEdgeData_m[iOpt], &selfConnData_m[iOpt], &vertexData_m[iOpt], paraData,
								&firstEdgeVertice[startIdx], &vertexNeighbor[startIdx],
								&data[iOpt]);
				}
			}
		} else
		{
		}
		athread_join();
	}
	for(i=0;i<BLOCKNUM64K;i++)
	{
		blockIdx = i*cpeBlockNum+BLOCKNUM64K;
		startIdx = blockStarts[4*blockIdx+2];
		blockIdx = (i+1)*cpeBlockNum-1;
		endIdx = blockStarts[4*blockIdx+3];

		for(int iOpt=0;iOpt<optNum;iOpt++)
		{
			backEdgeData_m[iOpt].fArraySizes = 0;
			frontEdgeData_m[iOpt].fArraySizes= endIdx-startIdx;
			selfConnData_m[iOpt].fArraySizes = 0;

			moveArrayPtr(frontEdgeData_m, frontEdgeArrays, startIdx);
			cOpt[iOpt].fun_host(&backEdgeData_m[iOpt], &frontEdgeData_m[iOpt], &selfConnData_m[iOpt], &vertexData_m[iOpt], paraData,
						&firstEdgeVertice[startIdx], &vertexNeighbor[startIdx],
						&data[iOpt]);
		}
	}

#undef moveArrayPtr
}

void MultiLevelBlockIterator::reorderEdgeDataUnsymm(Arrays* edgeData)
{
	for(int i=0;i<edgeData->fArrayNum;i++)
	{
		reorderEdgeArrayUnsymm(edgeData->floatArrays[i]);
	}
}

void MultiLevelBlockIterator::reorderEdgeArrayUnsymm(swFloat* array)
{
	swInt  edgeNumber    = this->getTopology()->getEdgeNumber();
	swInt  vertexNumber  = this->getTopology()->getVertexNumber();
	swInt  cpeBlockNum   = this->_cpeBlockNum;
	swInt  mshBlockNum   = this->_mshBlockNum;
	swInt* blockStarts  = this->getBlockStarts();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getTopology()->getStartVertices();
	swInt* neighbor     = this->getTopology()->getEndVertices();

	swFloat* tmp = (swFloat*)malloc(sizeof(swFloat)*edgeNumber*2);
	swInt* firstEdgeVertices = this->_firstEdgeVertices;
	swInt* vertexNeighbours  = this->_vertexNeighbours;
	swInt* accuVertexEdgeNumbers
		= this->getTopology()->getAccuVertexEdgeNumbers();
	swInt* accuStartVertexNumbers
		= this->getTopology()->getAccuStartVertexNumbers();
	swInt row,col,upIdx=0;
	swInt*  lowIdx = (swInt*)malloc(sizeof(swInt)*vertexNumber);
	swFloat* lower = (swFloat*)malloc(sizeof(swFloat)*edgeNumber);
	swFloat* upper = (swFloat*)malloc(sizeof(swFloat)*edgeNumber);
	for(int i=0;i<vertexNumber;i++) {lowIdx[i]=0;}
	for(int i=0;i<edgeNumber;i++) {lower[i]=-1;upper[i]=-1;}

	// COO->LDU
	for(int i=0;i<vertexNumber;i++)
	{
		for(int j=accuVertexEdgeNumbers[i];j<accuVertexEdgeNumbers[i+1];j++)
		{
			if(vertexNeighbours[j]>i)
			{
				upper[upIdx] = array[j];
				upIdx++;
			}else if(vertexNeighbours[j]<i)
			{
				assert(accuStartVertexNumbers[vertexNeighbours[j]]+lowIdx[vertexNeighbours[j]]<accuStartVertexNumbers[vertexNeighbours[j]+1]);
				lower[accuStartVertexNumbers[vertexNeighbours[j]]
					+lowIdx[vertexNeighbours[j]]] = array[j];
				lowIdx[vertexNeighbours[j]]++;
			}
		}
	}
//	for(int i=0;i<edgeNumber;i++)
//	{
//		if(upper[i]==-1 || lower[i]==-1)
//		{
//			printf("%d,%f,%f\n",i,lower[i],upper[i]);
//		}
//	}

	// LDU->reordered LDU
	Arrays backEdgeData, frontEdgeData;
	constructSingleArray(backEdgeData, 1, edgeNumber, COPYIN, lower);
	constructSingleArray(frontEdgeData, 1, edgeNumber, COPYIN, upper);
//	Arrays edgeData = {lower, upper, NULL, NULL, edgeNumber};
//	TODO
	reorderEdgeData(&backEdgeData, &frontEdgeData);
//	reorderEdgeData(&edgeData);

	// reordered LDU -> reordered COO data
	swInt* vertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* accuVertexEdgeNums
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);

	int blockIdx, startBlockIdx, endBlockIdx, row1, blockIdxV;
	swFloat value;
	int accuEdges= 0;
	int startIdx = 0;
	int endIdx   = 0;
	for(int i=0;i<mshBlockNum;i++)
	{
		for(int j=0;j<BLOCKNUM64K;j++)
		{
			int row = i*BLOCKNUM64K+j;
			int cellLen = vertexStarts[row+1]-vertexStarts[row];
			assert(cellLen<vertexNumber/cpeBlockNum*2);
			for(int col=0;col<cpeBlockNum;col++)
			{
				for(int k=0;k<vertexNumber/cpeBlockNum*2;k++)
				{
					vertexEdgeNumbers[k]=0;
					accuVertexEdgeNums[k]=0;
				}
				if(col==row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
		    			vertexEdgeNumbers[owner[k]-vertexStarts[row]]++;
		    			vertexEdgeNumbers[neighbor[k]-vertexStarts[row]]++;
		    		}
					for(int k=1;k<cellLen;k++)
					{
						accuVertexEdgeNums[k+1]
							= accuVertexEdgeNums[k]
							+ vertexEdgeNumbers[k-1];
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = neighbor[k]+1-vertexStarts[row];
//						if(accuVertexEdgeNums[row1]+startIdx==178894) printf("%d,%f\n",k,lower[k]);
						array[accuVertexEdgeNums[row1]+startIdx]=lower[k];
						accuVertexEdgeNums[row1]++;
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = owner[k]+1-vertexStarts[row];
//						if(accuVertexEdgeNums[row1]+startIdx==178894) printf("%d,%f\n",k,upper[k]);
						array[accuVertexEdgeNums[row1]+startIdx]=upper[k];
						accuVertexEdgeNums[row1]++;
					}
					startIdx += (endBlockIdx-startBlockIdx)*2;
				} else if(col < row)
				{
     				blockIdx = col*(1+2*this->_cpeBlockNum-col)/2+row-col;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
//						if(k-startBlockIdx+startIdx==178894) printf("%d,%f\n",k,lower[k]);
						array[k-startBlockIdx+startIdx] = lower[k];
		    		}
					startIdx += (endBlockIdx-startBlockIdx);
				} else if(col > row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2+col-row;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
//						if(k-startBlockIdx+startIdx==178894) printf("%d,%f\n",k,upper[k]);
						array[k-startBlockIdx+startIdx] = upper[k];
		    		}
					startIdx += (endBlockIdx-startBlockIdx);
				}
			}
		}
	}

	free(tmp);
	free(lower);
	free(upper);
	free(lowIdx);
	free(vertexEdgeNumbers);
	free(accuVertexEdgeNums);
}


void MultiLevelBlockIterator::reorderEdgeData(Arrays* backEdgeData,
			Arrays* frontEdgeData)
{
	map<swInt, swInt>::iterator iter;
	swFloat *tmpBED,*tmpFED;
	swInt edgeNum = frontEdgeData->fArraySizes;
	for(int i=0;i<backEdgeData->fArrayNum;i++)
	{
		swInt dims = backEdgeData->fArrayDims[i];
		tmpBED = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat));
		tmpFED = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat));
		for(iter=this->getEdgeMap().begin();
					iter!=this->getEdgeMap().end();iter++)
		{
//if(iter->second==1) printf("%d,%d,%d,%f,%f\n",iter->first,this->getTopology()->getStartVertices()[abs(iter->second)],this->getTopology()->getEndVertices()[abs(iter->second)],frontEdgeData->floatArrays[i][iter->first]),backEdgeData->floatArrays[i][iter->first];
			if(iter->second<0)
			{
				for(int j=0;j<dims;j++)
				{
					tmpBED[(-iter->second-1)*dims+j]
						= frontEdgeData->floatArrays[i][iter->first*dims+j];
					tmpFED[(-iter->second-1)*dims+j]
						= backEdgeData->floatArrays[i][iter->first*dims+j];
				}
			} else
			{
				for(int j=0;j<dims;j++)
				{
					tmpBED[(iter->second-1)*dims+j]
						= backEdgeData->floatArrays[i][iter->first*dims+j];
					tmpFED[(iter->second-1)*dims+j]
						= frontEdgeData->floatArrays[i][iter->first*dims+j];
				}
			}
		}
		for(int iedge=0;iedge<edgeNum;iedge++)
		{
			for(int j=0;j<dims;j++)
			{
				backEdgeData->floatArrays[i][iedge*dims+j]
					= tmpBED[iedge*dims+j];
				frontEdgeData->floatArrays[i][iedge*dims+j]
					= tmpFED[iedge*dims+j];
			}
		}
//printf("%f\n",frontEdgeData->floatArrays[0][0]);
		free(tmpBED);
		free(tmpFED);
	}
}

void MultiLevelBlockIterator::restoreEdgeData(Arrays* backEdgeData,
			Arrays* frontEdgeData)
{
	map<swInt, swInt>::iterator iter;
	swFloat *tmpBED,*tmpFED;
	swInt edgeNum = backEdgeData->fArraySizes;
	for(int i=0;i<backEdgeData->fArrayNum;i++)
	{
		swInt dims = backEdgeData->fArrayDims[i];
		tmpBED = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat));
		tmpFED = (swFloat*)malloc(edgeNum*dims*sizeof(swFloat));
		for(iter=this->getEdgeMap().begin();
					iter!=this->getEdgeMap().end();iter++)
		{
//if(iter->first==84) printf("%d,%f\n",iter->second,backEdgeData->floatArrays[i][(-iter->second-1)]);
//if(iter->first==56236) printf("%d,%d,%d,%f,%f\n",iter->second,this->getTopology()->getStartVertices()[abs(iter->second)],this->getTopology()->getEndVertices()[abs(iter->second)],frontEdgeData->floatArrays[i][iter->first]),backEdgeData->floatArrays[i][iter->first];
			if(iter->second<0)
			{
				for(int j=0;j<dims;j++)
				{
					tmpBED[iter->first*dims+j]
						= frontEdgeData->floatArrays[i][(-iter->second-1)*dims+j];
					tmpFED[iter->first*dims+j]
						= backEdgeData->floatArrays[i][(-iter->second-1)*dims+j];
				}
			} else
			{
				for(int j=0;j<dims;j++)
				{
					tmpBED[iter->first*dims+j]
						= backEdgeData->floatArrays[i][(iter->second-1)*dims+j];
					tmpFED[iter->first*dims+j]
						= frontEdgeData->floatArrays[i][(iter->second-1)*dims+j];
				}
			}
		}
		for(int iedge=0;iedge<edgeNum;iedge++)
		{
			for(int j=0;j<dims;j++)
			{
				backEdgeData->floatArrays[i][iedge*dims+j]
					= tmpBED[iedge*dims+j];
				frontEdgeData->floatArrays[i][iedge*dims+j]
					= tmpFED[iedge*dims+j];
			}
		}
		free(tmpBED);
		free(tmpFED);
	}
}

void MultiLevelBlockIterator::reorderVertexArray(swFloat* array)
{
	swFloat *tmp = (swFloat*)malloc
		(this->getTopology()->getVertexNumber()*(sizeof(swFloat)));
	map<swInt, swInt>::iterator iter;
	for(iter = this->getVertexMap().begin();
				iter!=this->getVertexMap().end();iter++)
	{
//		if(iter->first==0) cout<<iter->second<<endl;
		tmp[iter->second] = array[iter->first];
	}
	for(int i=0;i<this->getTopology()->getVertexNumber();i++)
	{
		array[i] = tmp [i];
	}
	free(tmp);
}

void MultiLevelBlockIterator::reorderVertexData(Arrays* vertexData)
{
	swFloat *tmp;
	swInt vertexNum = vertexData->fArraySizes;
	map<swInt, swInt>::iterator iter;
	for(int i=0;i<vertexData->fArrayNum;i++)
	{
		swInt dims = vertexData->fArrayDims[i];
		tmp = (swFloat*)malloc(vertexNum*dims*sizeof(swFloat));
		for(iter=this->getVertexMap().begin();
					iter!=this->getVertexMap().end();iter++)
		{
			for(int j=0;j<dims;j++)
			{
				tmp[iter->second*dims+j]
					= vertexData->floatArrays[i][iter->first*dims+j];
			}
		}
		for(int ivertex=0;ivertex<vertexNum;ivertex++)
		{
			for(int j=0;j<dims;j++)
			{
				vertexData->floatArrays[i][ivertex*dims+j]
					= tmp[ivertex*dims+j];
			}
		}
		free(tmp);
	}
}

void MultiLevelBlockIterator::initOwnNeiSendList()
{
	swInt mshBlockNum = this->_mshBlockNum;
	swInt cpeBlockNum = this->_cpeBlockNum;
	swInt *owner,*neighbor,*cellStarts,*blockStarts,*ownNeiSendList;
	owner = this->getTopology()->getStartVertices();
	neighbor = this->getTopology()->getEndVertices();
	cellStarts = this->_vertexStarts;
	blockStarts = this->_blockStarts;

	vector<vector<int> > dataList(cpeBlockNum);
	ownNeiSendList
		= (swInt*)malloc(cpeBlockNum*BLOCKNUM64K*sizeof(swInt));
	for(int i=0;i<cpeBlockNum*BLOCKNUM64K;i++)
	{
		ownNeiSendList[i] = 0;
	}

	// Get the count of edges of each column block
	for(int irow=0;irow<cpeBlockNum;irow++)
	{
		int iseg = irow/BLOCKNUM64K;
		int localRow = irow-iseg*BLOCKNUM64K;
		for(int icol=localRow;icol<BLOCKNUM64K;icol++)
		{
			int idx = irow*BLOCKNUM64K+icol;
			int symmIdx = (icol+iseg*BLOCKNUM64K)*BLOCKNUM64K+localRow;
			int blockIdx
				= irow*(1+2*cpeBlockNum-irow)/2+icol+iseg*BLOCKNUM64K-irow;
			ownNeiSendList[idx]
				= blockStarts[4*blockIdx+3]-blockStarts[4*blockIdx+2];
			ownNeiSendList[symmIdx] = ownNeiSendList[idx];
		}
	}

	for(int irow=0;irow<cpeBlockNum;irow++)
	{
		int iseg = irow/BLOCKNUM64K;
		int localRow = irow-iseg*BLOCKNUM64K;
		for(int icol=0;icol<BLOCKNUM64K;icol++)
		{
			int idx = irow*BLOCKNUM64K+icol;
			if(icol==localRow)
			{
				ownNeiSendList[idx] = 0;
				continue;
			}
			if(ownNeiSendList[idx]==0) continue;
			int edgeNum = ownNeiSendList[idx];
			while(edgeNum>6)
			{
				dataList[irow].push_back(6);
				edgeNum -= 6;
			}
			dataList[irow].push_back(edgeNum);
			ownNeiSendList[idx] = (ownNeiSendList[idx]-1)/6+1;
		}
	}

	//printf("start Init sendList...\n");
	this->_schedule_data = new RlmpiInfo[mshBlockNum];
	for(int iseg=0;iseg<mshBlockNum;iseg++)
	{
		RlmpiInitializer reg;
		vector<vector<int> > sendList(BLOCKNUM64K);
		vector<vector<int> > nDataList(BLOCKNUM64K);
		for(int localRow=0;localRow<BLOCKNUM64K;localRow++)
		{
			int irow = iseg*BLOCKNUM64K+localRow;
			int packIdx = 0;
			for(int icol=0;icol<BLOCKNUM64K;icol++)
			{
				int idx = irow*BLOCKNUM64K+icol;
//printf("iseg: %d, sendId: %d, recvId: %d, sendNum: %d\n",iseg,localRow,icol,ownNeiSendList[idx]);
				for(int k=0;k<ownNeiSendList[idx];k++)
				{
//if(localRow==9 && iseg==9) printf("%d,%d,%d,%d,%d\n",irow,icol,packIdx,dataList[irow][packIdx],ownNeiSendList[idx]);
					sendList[localRow].push_back(icol);
					nDataList[localRow].push_back(dataList[irow][packIdx]);
					packIdx++;
				}
			}
		}
		reg.init(sendList, nDataList);
		reg.copyRlmpiInfo(&this->_schedule_data[iseg]);
	}
	//printf("finish Init sendList!\n");
}

void MultiLevelBlockIterator::writeTopology()
{
	swInt edgeNum   =this->getTopology()->getEdgeNumber();
	swInt *owner    = this->getTopology()->getStartVertices();
	swInt *neighbor = this->getTopology()->getEndVertices();
	printf("start writing file...\n");
	std::fstream file("owner_MLB",ios::out);
	if(!file) return;
	for(int i=0;i<edgeNum;i++)
	{
		printf("%d\n",i);
		file<<owner[i]<<std::endl;
	}
	file.close();
	printf("start writing file...\n");

	file.open("neighbour_MLB",ios::out);
	if(!file) return;
	for(int i=0;i<edgeNum;i++)
	{
		file<<neighbor[i]<<std::endl;
	}
	file.close();
}


}

