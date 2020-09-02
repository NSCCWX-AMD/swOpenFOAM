#include "rowSubsectionIterator.hpp"
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <sys/time.h>
#include <fstream>
#include "rowSubsectionIterator.h"
#include "rowSubsections.hpp"
#include "rowSubsection.h"

#ifdef __cplusplus
extern "C"{
#endif
#include <athread.h>
void slave_rowSubsectionIterator_v2e_slave(RSS_vertex2EdgePara *para);
void slave_rowSubsectionIterator_e2v_slave(RSS_vertex2EdgePara *para);
void slave_rowSubsectionIterator_array_slave(RSS_vertex2EdgePara *para);
//void slave_multiLevelBlockIterator_v2e_slave(MLB_vertex2EdgePara *para);
//void slave_multiLevelBlockIterator_array_slave(MLB_arrayPara *para);
//#include "edge2VertexIter_host.h"
//#include "vertex2EdgeIter_host.h"
#ifdef __cplusplus
}
#endif

#define LDM_MIN (64*1024*0.90)
#define LDM_MAX (64*1024*0.94)

namespace UNAT
{

// ***********************Constructor*************************************
RowSubsectionIterator::RowSubsectionIterator(coupledOperator* opt, 
			Topology &topo,	swInt* cellWeights, swInt* edgeWeights, 
			int pattern, bool duplicate)
		: Iterator(topo, cellWeights, edgeWeights, duplicate)
{
	swInt fSize, cSize, *lPtr, *uPtr;
	if(pattern==V2E)
	{
		fSize = this->getTopology()->getEdgeNumber()*2;
		cSize = this->getTopology()->getVertexNumber();
		lPtr = this->getTopology()->getFirstEdgeVertices();
		uPtr = this->getTopology()->getVertexNeighbours();
	} else if(pattern==E2V)
	{
		fSize = this->getTopology()->getEdgeNumber();
		cSize = this->getTopology()->getVertexNumber();
		lPtr = this->getTopology()->getStartVertices();
		uPtr = this->getTopology()->getEndVertices();
	} else if(pattern==ARRAY)
	{
		return;
	} else
	{
		printf("\nError: pattern is not supported!\n");
	}

	swInt maxFaces = (swInt)(64*1024/
				(sizeof(swInt)*2+edgeWeights[0]*sizeof(swFloat)));
	this->maxFacesInSeg_ = (swInt*)malloc(BLOCKNUM64K*sizeof(swInt));
	this->maxCellsInSeg_ = (swInt*)malloc(BLOCKNUM64K*sizeof(swInt));
	this->maxColsInSeg_ = (swInt*)malloc(BLOCKNUM64K*sizeof(swInt));
	this->maxColNum_ = 0;

	// Get datasize for row-subsections
	// choose the maximum dimension of vertexData with COPYOUT attribute.
	Arrays* vertexData = opt->data->vertexData;
	int dataSize = 0;
	for(int iArray=0;iArray<vertexData->fArrayNum;iArray++)
	{
		if(vertexData->fArrayInOut[iArray]==COPYOUT)
		{
			dataSize = dataSize > vertexData->fArrayDims[iArray] ?
				dataSize : vertexData->fArrayDims[iArray];
		}
	}

	float coe = 1.0;
	int istep = 0;
	int totalTotalCols = 0;
	while(true)
	{
		istep++;
		int maxTotalLength = 0;
		maxFaces = (swInt)(maxFaces*coe);
		this->secs_ = new SWFoam::RowSubsections(fSize, cSize, BLOCKNUM64K,
					4, lPtr, uPtr, maxFaces, sizeof(swFloat));
		const rowSubsection** secs = this->secs_->getSubsections();
		swInt secNumInSeg = this->secs_->getSecNumInSeg();
		totalTotalCols = 0;
		maxColNum_ = 0;
		for(int iseg=0;iseg<BLOCKNUM64K;iseg++)
		{
			this->maxFacesInSeg_[iseg] = 0;
			this->maxCellsInSeg_[iseg] = 0;
			this->maxColsInSeg_[iseg] = 0;
			for(int isec=0;isec<secNumInSeg;isec++)
			{
				const rowSubsection subsection  = secs[iseg][isec];
				const swInt faceStart           = subsection.faceStart;
				const swInt nFaces              = subsection.nFaces;
				const swInt cellStart           = subsection.cellStart;
				const swInt nCells              = subsection.nCells;
				const swInt nSecs               = subsection.nSecs;
				const swInt *colStartsAndCounts
					= subsection.colStartsAndCounts;
				swInt totalCols = 0;
				swInt totalColIdx = nSecs;
				while(totalColIdx!=0)
				{
					totalCols += colStartsAndCounts[2*totalColIdx-1];
					maxColNum_ =
						maxColNum_ > colStartsAndCounts[2*totalColIdx-1] ?
						maxColNum_ : colStartsAndCounts[2*totalColIdx-1];
					totalColIdx--;
				}
				totalTotalCols += totalCols;
				maxColsInSeg_[iseg] = maxColsInSeg_[iseg] > totalCols
					? maxColsInSeg_[iseg] : totalCols;
				maxFacesInSeg_[iseg] = maxFacesInSeg_[iseg] > nFaces
					? maxFacesInSeg_[iseg] : nFaces;
//				swInt nCells = lPtr[faceStart+nFaces-1]-lPtr[faceStart]+1;
				maxCellsInSeg_[iseg] = maxCellsInSeg_[iseg] > nCells
					? maxCellsInSeg_[iseg] : nCells;
			}
			int totalLength
				// owner, neighbor
				= maxFacesInSeg_[iseg]*2*sizeof(swInt)
				// Faces
				+ maxFacesInSeg_[iseg]*edgeWeights[0]*sizeof(swFloat)
				// Cells
				+ (maxColsInSeg_[iseg]+maxCellsInSeg_[iseg])
				* (cellWeights[0])*sizeof(swFloat)
				+ 1024*dataSize;
			maxTotalLength = maxTotalLength > totalLength
				? maxTotalLength : totalLength;
		}
		// printf("istep: %d, coe: %f, LDM_size: %d\n", istep, coe, maxTotalLength);
		if(maxTotalLength <= LDM_MAX && maxTotalLength >= LDM_MIN)
		{
			break;
		} else if(maxTotalLength <LDM_MIN)
		{
			if(istep>5) break;
			coe = (double)(LDM_MIN/maxTotalLength);
		} else
		{
			coe = (double)(LDM_MAX/maxTotalLength);
		}
		delete(this->secs_);
	}
	// printf("RowSubsection--secInSeg: %d\n",this->secs_->getSecNumInSeg());
	// printf("nFaces:%d, nCells:%d, totalCols: %d, DMA_size:%f GB\n",fSize, cSize, totalTotalCols,(double)((fSize*edgeWeights[0]+cSize*(cellWeights[0]+1)+totalTotalCols*cellWeights[0]/2)*sizeof(swFloat)+fSize*2*sizeof(swInt))/1024/1024/1024);
}



void RowSubsectionIterator::vertex2EdgeIteration(Arrays *paraData,
			coupledOperator *cOpt, swInt optNum)
{
	swInt *firstEdgeVertices = this->getTopology()->getFirstEdgeVertices();
	swInt *vertexNeighbours  = this->getTopology()->getVertexNeighbours();
	swInt cellNum = this->getTopology()->getVertexNumber();
	swInt edgeNum = this->getTopology()->getEdgeNumber();

printf("colRoundNum: %d\n",this->secs_->getColRoundNum());
	RSS_vertex2EdgePara para = {
		firstEdgeVertices,
		vertexNeighbours,
		this->secs_->getSubsections(),
		this->secs_->getSecNumInSeg(),
		this->secs_->getColRoundNum(),
		cellNum,
		edgeNum,
		this->maxFacesInSeg_,
		this->maxCellsInSeg_,
		this->maxColsInSeg_,
		this->maxColNum_,

		paraData,

		cOpt,
		optNum
	};
	__real_athread_spawn((void*)slave_rowSubsectionIterator_v2e_slave,
				&para);
	athread_join();
}

void RowSubsectionIterator::edge2VertexIteration(Arrays *paraData,
			coupledOperator *cOpt, swInt optNum)
{
	swInt *startVertices = this->getTopology()->getStartVertices();
	swInt *endVertices   = this->getTopology()->getEndVertices();
	swInt cellNum = this->getTopology()->getVertexNumber();
	swInt edgeNum = this->getTopology()->getEdgeNumber();

	RSS_vertex2EdgePara para = {
		startVertices,
		endVertices,
		this->secs_->getSubsections(),
		this->secs_->getSecNumInSeg(),
		this->secs_->getColRoundNum(),
		cellNum,
		edgeNum,
		this->maxFacesInSeg_,
		this->maxCellsInSeg_,
		this->maxColsInSeg_,
		this->maxColNum_,

		paraData,

		cOpt,
		optNum
	};

	__real_athread_spawn((void*)slave_rowSubsectionIterator_e2v_slave,
				&para);
	athread_join();
}

void RowSubsectionIterator::arrayIteration(Arrays *paraData,
			coupledOperator *cOpt, swInt optNum)
{
	swInt edgeNum = cOpt->data->frontEdgeData[0].fArraySizes;

	RSS_vertex2EdgePara para = {
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		0,
		edgeNum,
		NULL,
		NULL,
		NULL,
		NULL,

		paraData,

		cOpt,
		optNum
	};

	__real_athread_spawn((void*)slave_rowSubsectionIterator_array_slave,
				&para);
	athread_join();
}

}
