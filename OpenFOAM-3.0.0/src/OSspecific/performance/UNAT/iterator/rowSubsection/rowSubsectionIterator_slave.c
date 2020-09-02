#include "swMacro.h"
#include "slave.h"
#include <assert.h>
#include "slaveUtils.h"
#include "rowSubsectionIterator.h"


#define copyArrayWithoutData(srcArray, dstArray, size) \
{ \
	int iArray; \
	char* tmp_ptr = ldm_space_end; \
	ldm_space_end = ldm_ptr; \
	(dstArray).fArraySizes = (size); \
	(dstArray).fArrayNum   = (srcArray)->fArrayNum; \
	if((dstArray).fArrayNum>0) \
	{ \
		LDM_NEW((dstArray).fArrayInOut,swInt,(dstArray).fArrayNum); \
		LDM_NEW((dstArray).fArrayDims, swInt,(dstArray).fArrayNum); \
		LDM_NEW((dstArray).floatArrays,swFloat*,(dstArray).fArrayNum); \
		DMA_Get((dstArray).fArrayInOut, (srcArray)->fArrayInOut, \
					(dstArray).fArrayNum*sizeof(swInt)); \
		DMA_Get((dstArray).fArrayDims, (srcArray)->fArrayDims, \
					(dstArray).fArrayNum*sizeof(swInt)); \
	} \
	ldm_ptr = ldm_space_end; \
	ldm_space_end = ldm_space_end>tmp_ptr?ldm_space_end:tmp_ptr; \
}

#define copyArrayWithoutData_new(srcArray, dstArray, size) \
{ \
	int iArray; \
	(dstArray).fArraySizes = (size); \
	(dstArray).fArrayNum   = (srcArray)->fArrayNum; \
	if((dstArray).fArrayNum>0) \
	{ \
		LDM_NEW((dstArray).fArrayInOut,swInt,(dstArray).fArrayNum); \
		LDM_NEW((dstArray).fArrayDims, swInt,(dstArray).fArrayNum); \
		LDM_NEW((dstArray).floatArrays,swFloat*,(dstArray).fArrayNum); \
		DMA_Get((dstArray).fArrayInOut, (srcArray)->fArrayInOut, \
					(dstArray).fArrayNum*sizeof(swInt)); \
		DMA_Get((dstArray).fArrayDims, (srcArray)->fArrayDims, \
					(dstArray).fArrayNum*sizeof(swInt)); \
	} \
}


void rowSubsectionIterator_v2e_slave(RSS_vertex2EdgePara *para)
{
	INIT_LDM_SPACE(62000);
	RSS_vertex2EdgePara para_s = *para;
	const swInt *lPtr          = para_s.firstEdgeVertice;
	const swInt *uPtr          = para_s.vertexNeighbor;
	const swInt faceNum        = para_s.edgeNumber;
	const swInt cellNum        = para_s.vertexNumber;
	const swInt secNumInSeg    = para_s.secNumInSeg;
	const struct rowSubsection** secs = para_s.secs;
	const swInt maxFaces       = para_s.maxFacesInSeg[_MYID];
	const swInt maxCells       = para_s.maxCellsInSeg[_MYID];
	const swInt maxCols        = para_s.maxColsInSeg[_MYID];
	
	swInt *owner, *neighbor;
	LDM_NEW(owner,    swInt, maxFaces);
	LDM_NEW(neighbor, swInt, maxFaces);

	int isec, iOpt;
	int optNum = para_s.optNum;
	Arrays frontEdgeData_m[optNum], backEdgeData_m[optNum];
   	Arrays selfConnData_m[optNum], vertexData_m[optNum];
	Arrays frontEdgeData[optNum],   backEdgeData[optNum];
	Arrays selfConnData[optNum],   vertexData[optNum];
	FieldData data_s[optNum];
	for(iOpt=0;iOpt<optNum;iOpt++)
	{
		data_s[iOpt].frontEdgeData = &frontEdgeData[iOpt];
		data_s[iOpt].backEdgeData  = &backEdgeData[iOpt];
		data_s[iOpt].selfConnData  = &selfConnData[iOpt];
		data_s[iOpt].vertexData    = &vertexData[iOpt];
	}
	swFloat *edgePtr[10], *vertexPtr[10], *selfConnPtr[10];
	int eMax=0,vMax=0,sMax=0;
	//char *ldm_tag = ldm_space_end;
	//ldm_space_end += 1024;
	for(isec=0;isec<secNumInSeg;isec++)
	{
		const struct rowSubsection subsection = secs[_MYID][isec];
		const swInt faceStart          = subsection.faceStart;
		const swInt nFaces             = subsection.nFaces;
		const swInt nSecs              = subsection.nSecs;
		const swInt colSAndC[nSecs*2];

		DMA_Status status=0;
		DMA_IGet(colSAndC, subsection.colStartsAndCounts,
					nSecs*2*sizeof(swInt), &status);
		DMA_IGet(owner,    &lPtr[faceStart], nFaces*sizeof(swInt), &status);
		DMA_IGet(neighbor, &uPtr[faceStart], nFaces*sizeof(swInt), &status);
		DMA_Wait(&status, 3);
		const swInt cellStart = owner[0];
		const swInt nCells    = owner[nFaces-1] - owner[0] + 1;

		int iArray;
		for(iOpt=0;iOpt<optNum;iOpt++)
		{
			int eIdx=0,vIdx=0,sIdx=0;

			FieldData data = *para_s.cOpt[iOpt].data;
			e2v_slaveFunPtr fun_slave = para_s.cOpt[iOpt].fun_slave;

			if(isec==0)
			{
				// Get arrays data from MPE but floatArrays
				status = 0;
				DMA_IGet(&frontEdgeData_m[iOpt], data.frontEdgeData, 
							sizeof(Arrays),	&status);
				DMA_IGet(&vertexData_m[iOpt], data.vertexData, sizeof(Arrays),
							&status);
				DMA_IGet(&selfConnData_m[iOpt], data.selfConnData, 
							sizeof(Arrays), &status);
				DMA_Wait(&status, 3);
				//ldm_space_end = ldm_tag;
				copyArrayWithoutData_new(&frontEdgeData_m[iOpt], 
							frontEdgeData[iOpt], nFaces);
				copyArrayWithoutData_new(&vertexData_m[iOpt], vertexData[iOpt], 
							nCells);
				copyArrayWithoutData_new(&selfConnData_m[iOpt], 
							selfConnData[iOpt], nCells);
				// backEdgeData
				constructEmptyArray(backEdgeData[iOpt]);
			}

			status = 0;
			int dma_num = 0;

			// FrontEdgeData
			for(iArray=0;iArray<frontEdgeData[iOpt].fArrayNum;iArray++)
			{
				int dims = frontEdgeData[iOpt].fArrayDims[iArray];
				if(frontEdgeData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					frontEdgeData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
					continue;
				} else if(eIdx<eMax)
				{
					frontEdgeData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
				} else
				{
					LDM_NEW(frontEdgeData[iOpt].floatArrays[iArray], 
								swFloat, maxFaces*dims);
					edgePtr[eIdx++] 
						= frontEdgeData[iOpt].floatArrays[iArray];
					eMax++;
				}
				DMA_IGet(frontEdgeData[iOpt].floatArrays[iArray],
							&frontEdgeData_m[iOpt].floatArrays
							[iArray][faceStart*dims],
							nFaces*dims*sizeof(swFloat),&status);
				dma_num++;
			}


			// vertexData
			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
			{
				int dims = vertexData[iOpt].fArrayDims[iArray];
				if(vertexData[iOpt].fArrayInOut[iArray]==COPYIN) continue;
				if(vertexData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					vertexData[iOpt].floatArrays[iArray]
						= vertexPtr[vIdx++];
					continue;
				} else if(vIdx<vMax)
				{
					vertexData[iOpt].floatArrays[iArray] 
						= vertexPtr[vIdx++];
				} else
				{
					LDM_NEW(vertexData[iOpt].floatArrays[iArray], swFloat,
								maxCells*dims);
					vertexPtr[vIdx++] 
						= vertexData[iOpt].floatArrays[iArray];
					vMax++;
				}
				DMA_IGet(vertexData[iOpt].floatArrays[iArray],
							&vertexData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							nCells*dims*sizeof(swFloat),
							&status);
				dma_num++;
			}
			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
			{
				int dims = vertexData[iOpt].fArrayDims[iArray];
				if(vertexData[iOpt].fArrayInOut[iArray]==COPYOUT) continue;
				if(vertexData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					vertexData[iOpt].floatArrays[iArray] 
						= vertexPtr[vIdx++];
					continue;
				} else if(vIdx<vMax)
				{
					vertexData[iOpt].floatArrays[iArray] 
						= vertexPtr[vIdx++];
				} else
				{
					LDM_NEW(vertexData[iOpt].floatArrays[iArray], swFloat,
								(maxCols+maxCells)*dims);
					vertexPtr[vIdx++] 
						= vertexData[iOpt].floatArrays[iArray];
					vMax++;
				}
				int tmpStart = maxCells;
				int icol;
				for(icol=0;icol<nSecs;icol++)
				{
					swInt colStart = colSAndC[2*icol];
					swInt colNum   = colSAndC[2*icol+1];
					DMA_IGet(&vertexData[iOpt].floatArrays
							[iArray][tmpStart*dims],
							&vertexData_m[iOpt].floatArrays
							[iArray][colStart*dims],
							colNum*dims*sizeof(swFloat), &status);
					dma_num++;
					tmpStart += colNum;
				}
				DMA_IGet(vertexData[iOpt].floatArrays[iArray],
							&vertexData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							nCells*dims*sizeof(swFloat), &status);
				dma_num++;
//if(_MYID==38) printf("%p,%f,%f\n",ldm_space_end,vertexData.floatArrays[iArray][maxCols+256],vertexData_m.floatArrays[iArray][30120]);
			}

			// selfConnData
			for(iArray=0;iArray<selfConnData[iOpt].fArrayNum;iArray++)
			{
				int dims = selfConnData[iOpt].fArrayDims[iArray];
				if(selfConnData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					selfConnData[iOpt].floatArrays[iArray] 
						= selfConnPtr[sIdx++];
					continue;
				} else if(sIdx<sMax)
				{
					selfConnData[iOpt].floatArrays[iArray] 
						= selfConnPtr[sIdx++];
				} else if(eIdx<eMax)
				{
					selfConnData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
				} else
				{
					LDM_NEW(selfConnData[iOpt].floatArrays[iArray], swFloat,
								maxCells*dims);
					selfConnPtr[sIdx++] 
						= selfConnData[iOpt].floatArrays[iArray];
					sMax++;
				}
				DMA_IGet(selfConnData[iOpt].floatArrays[iArray],
							&selfConnData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							nCells*dims*sizeof(swFloat), &status);
				dma_num++;
			}

			// Get neighbor
			int iface,icell;
			for(iface=0;iface<nFaces;iface++)
			{
				owner[iface] -= cellStart;
				swInt col_global = neighbor[iface];
				swInt col_local  = 0;
				int isec1=0;
				while(col_global<colSAndC[2*isec1] ||
							col_global>=colSAndC[2*isec1+1]+colSAndC[2*isec1])
				{
					if(isec1<nSecs)
					{
						col_local += colSAndC[2*isec1+1];
						isec1++;
					} else
					{
						printf("Error at new neighbor!\n");
					}
				}
				neighbor[iface] = col_local+col_global-colSAndC[2*isec1]
					+ maxCells;
			}

			DMA_Wait(&status, dma_num);

			// Compute

			frontEdgeData[iOpt].fArraySizes = 0;
			selfConnData[iOpt].fArraySizes  = nCells;
			fun_slave(NULL,NULL,NULL,NULL,NULL,
						owner,neighbor,&data_s[iOpt]);
//			swFloat *x = accessArray(&vertexData, 0);
//			swFloat *b = accessArray(&vertexData, 1);
//			swFloat *edge = accessArray(&frontEdgeData, 0);
//			swFloat *diag = accessArray(&selfConnData, 0);
//
//				for(icell=0;icell<nCells;icell++)
//				{
////if(icell+cellStart==30120) printf("%d,%d,%d,%f,%f,%f\n",_MYID,icell,cellStart,b[icell],diag[icell],x[icell]);
//					b[icell] += diag[icell]*x[icell];
//				}

			frontEdgeData[iOpt].fArraySizes = nFaces;
			selfConnData[iOpt].fArraySizes  = 0;
			fun_slave(NULL,NULL,NULL,NULL,NULL,
						owner,neighbor,&data_s[iOpt]);

//			x = accessArray(&vertexData, 0);
//			b = accessArray(&vertexData, 1);
//			for(iface=0;iface<nFaces;iface++)
//			{
//				b[owner[iface]] += edge[iface]*x[neighbor[iface]];
//			}

			// vertexData(output)
			dma_num = 0;
			status = 0;
			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
			{
				int dims = vertexData[iOpt].fArrayDims[iArray];
				if(vertexData[iOpt].fArrayInOut[iArray]==COPYIN) continue;
				DMA_IPut(&vertexData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							vertexData[iOpt].floatArrays[iArray],
							nCells*dims*sizeof(swFloat),&status);
				dma_num++;
			}
//			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
//			{
//				int dims = vertexData[iOpt].fArrayDims[iArray];
//				if(vertexData[iOpt].fArrayInOut[iArray]==COPYIN) continue;
//				int tmpStart = 0;
//				int icol;
//				for(icol=0;icol<nSecs;icol++)
//				{
//					swInt colStart = colSAndC[2*icol];
//					swInt colNum   = colSAndC[2*icol+1];
//					DMA_IPut(&vertexData_m[iOpt].floatArrays
//							[iArray][colStart*dims],
//							&vertexData[iOpt].floatArrays
//							[iArray][tmpStart*dims],
//							colNum*dims*sizeof(swFloat), &status);
//					dma_num++;
//					tmpStart += colNum;
//				}
//			}
			DMA_Wait(&status, dma_num);
		}
	}
}

void rowSubsectionIterator_e2v_slave(RSS_vertex2EdgePara *para)
{
	INIT_LDM_SPACE(61500);
	RSS_vertex2EdgePara para_s = *para;
	const swInt *lPtr          = para_s.firstEdgeVertice;
	const swInt *uPtr          = para_s.vertexNeighbor;
	const swInt faceNum        = para_s.edgeNumber;
	const swInt cellNum        = para_s.vertexNumber;
	const swInt secNumInSeg    = para_s.secNumInSeg;
	const struct rowSubsection** secs = para_s.secs;
	const swInt maxFaces       = para_s.maxFacesInSeg[_MYID];
	const swInt maxCells       = para_s.maxCellsInSeg[_MYID];
	const swInt maxCols        = para_s.maxColsInSeg[_MYID];
	const swInt colRoundNum    = para_s.colRoundNum;
	const swInt maxColNum      = para_s.maxColNum;
		
	int iArray;	
	Arrays paraData,paraData_m;
	DMA_Get(&paraData_m, para_s.paraData, sizeof(Arrays));
	copyArrayWithoutData_new(&paraData_m, paraData, paraData_m.fArraySizes);
	for(iArray=0;iArray<paraData.fArrayNum;iArray++)
	{
		int length = paraData.fArraySizes;
		LDM_NEW(paraData.floatArrays[iArray], swFloat, length);
		DMA_Get(paraData.floatArrays[iArray],
			paraData_m.floatArrays[iArray],
			length*sizeof(swFloat));
	}

	swInt *owner, *neighbor;
	LDM_NEW(owner,    swInt, maxFaces);
	LDM_NEW(neighbor, swInt, maxFaces);

	int isec, iOpt;
	int optNum = para_s.optNum;
	Arrays frontEdgeData_m[optNum], backEdgeData_m[optNum];
   	Arrays selfConnData_m[optNum], vertexData_m[optNum];
	Arrays frontEdgeData[optNum],   backEdgeData[optNum];
	Arrays selfConnData[optNum],   vertexData[optNum];
	FieldData data_s[optNum];
	for(iOpt=0;iOpt<optNum;iOpt++)
	{
		data_s[iOpt].frontEdgeData = &frontEdgeData[iOpt];
		data_s[iOpt].backEdgeData  = &backEdgeData[iOpt];
		data_s[iOpt].selfConnData  = &selfConnData[iOpt];
		data_s[iOpt].vertexData    = &vertexData[iOpt];
	}
	swFloat *edgePtr[10], *vertexPtr[10], *selfConnPtr[10];
	int eMax=0,vMax=0,sMax=0;
	int maxDims=0;
	swFloat *tmpData;
//	swInt cellStart, nCells;
	//char *ldm_tag = ldm_space_end;
	//ldm_space_end += 1024;
	for(isec=0;isec<secNumInSeg;isec++)
	{
		const struct rowSubsection subsection = secs[_MYID][isec];
		const swInt faceStart          = subsection.faceStart;
		const swInt nFaces             = subsection.nFaces;
		const swInt cellStart          = subsection.cellStart;
		const swInt nCells             = subsection.nCells;
		const swInt nSecs              = subsection.nSecs;
		const swInt colRound           = subsection.colRound;
		const swInt colSAndC[nSecs*2];


		DMA_Status status=0;
		DMA_IGet(colSAndC, subsection.colStartsAndCounts,
					nSecs*2*sizeof(swInt), &status);
		DMA_IGet(owner,    &lPtr[faceStart], nFaces*sizeof(swInt), &status);
		DMA_IGet(neighbor, &uPtr[faceStart], nFaces*sizeof(swInt), &status);
		DMA_Wait(&status, 3);
// if(_MYID==0) printf("(faceStart:) %d\n", faceStart,lPtr[faceStart],owner[0]);
//		if(isec==0) cellStart = owner[0];
//		else cellStart += nCells;
//		nCells    = owner[nFaces-1] - cellStart + 1;


		for(iOpt=0;iOpt<optNum;iOpt++)
		{
			int eIdx=0,vIdx=0,sIdx=0;

			FieldData data = *para_s.cOpt[iOpt].data;
			e2v_slaveFunPtr fun_slave = para_s.cOpt[iOpt].fun_slave;

			if(isec==0)
			{
				// Get arrays data from MPE but floatArrays
				status = 0;
				DMA_IGet(&frontEdgeData_m[iOpt], data.frontEdgeData, 
							sizeof(Arrays),	&status);
				DMA_IGet(&backEdgeData_m[iOpt], data.backEdgeData, 
							sizeof(Arrays),	&status);
				DMA_IGet(&vertexData_m[iOpt], data.vertexData,
							sizeof(Arrays), &status);
				DMA_IGet(&selfConnData_m[iOpt], data.selfConnData, 
							sizeof(Arrays), &status);
				DMA_Wait(&status, 4);
				//ldm_space_end = ldm_tag;
				copyArrayWithoutData_new(&frontEdgeData_m[iOpt], 
							frontEdgeData[iOpt], nFaces);
				copyArrayWithoutData_new(&backEdgeData_m[iOpt], 
							backEdgeData[iOpt], nFaces);
				copyArrayWithoutData_new(&vertexData_m[iOpt],
							vertexData[iOpt], nCells);
				copyArrayWithoutData_new(&selfConnData_m[iOpt], 
							selfConnData[iOpt], nCells);
			}

			status = 0;
			int dma_num = 0;

			// FrontEdgeData
			for(iArray=0;iArray<frontEdgeData[iOpt].fArrayNum;iArray++)
			{
				int dims = frontEdgeData[iOpt].fArrayDims[iArray];
				if(frontEdgeData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					frontEdgeData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
					continue;
				} else if(eIdx<eMax)
				{
					frontEdgeData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
				} else
				{
					LDM_NEW(frontEdgeData[iOpt].floatArrays[iArray], 
								swFloat, maxFaces*dims);
					edgePtr[eIdx++] 
						= frontEdgeData[iOpt].floatArrays[iArray];
					eMax++;
				}
				DMA_IGet(frontEdgeData[iOpt].floatArrays[iArray],
							&frontEdgeData_m[iOpt].floatArrays
							[iArray][faceStart*dims],
							nFaces*dims*sizeof(swFloat),&status);
				dma_num++;
			}

			// BackEdgeData
			for(iArray=0;iArray<backEdgeData[iOpt].fArrayNum;iArray++)
			{
				int dims = backEdgeData[iOpt].fArrayDims[iArray];
				if(backEdgeData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					backEdgeData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
					continue;
				} else if(eIdx<eMax)
				{
					backEdgeData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
				} else
				{
					LDM_NEW(backEdgeData[iOpt].floatArrays[iArray], 
								swFloat, maxFaces*dims);
					edgePtr[eIdx++] 
						= backEdgeData[iOpt].floatArrays[iArray];
					eMax++;
				}
				DMA_IGet(backEdgeData[iOpt].floatArrays[iArray],
							&backEdgeData_m[iOpt].floatArrays
							[iArray][faceStart*dims],
							nFaces*dims*sizeof(swFloat),&status);
				dma_num++;
			}


			// vertexData
			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
			{
				int dims = vertexData[iOpt].fArrayDims[iArray];
				if(vertexData[iOpt].fArrayInOut[iArray]!=COPYOUT) continue;
				if(vertexData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					vertexData[iOpt].floatArrays[iArray]
						= vertexPtr[vIdx++];
					continue;
				} else if(vIdx<vMax)
				{
					vertexData[iOpt].floatArrays[iArray] 
						= vertexPtr[vIdx++];
				} else
				{
					LDM_NEW(vertexData[iOpt].floatArrays[iArray], swFloat,
								(maxCells+maxCols)*dims);
					vertexPtr[vIdx++] 
						= vertexData[iOpt].floatArrays[iArray];
					vMax++;
				}
				DMA_IGet(vertexData[iOpt].floatArrays[iArray],
							&vertexData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							nCells*dims*sizeof(swFloat),
							&status);
				dma_num++;
				maxDims = MAX(maxDims, dims);
				int icol;
				for(icol=maxCells*dims;icol<(maxCells+maxCols)*dims;icol++)
				{
					vertexData[iOpt].floatArrays[iArray][icol] = 0;
				}
			}
			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
			{
				int dims = vertexData[iOpt].fArrayDims[iArray];
				if(vertexData[iOpt].fArrayInOut[iArray]==COPYOUT) continue;
				if(vertexData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					vertexData[iOpt].floatArrays[iArray] 
						= vertexPtr[vIdx++];
					continue;
				} else if(vIdx<vMax)
				{
					vertexData[iOpt].floatArrays[iArray] 
						= vertexPtr[vIdx++];
				} else
				{
					LDM_NEW(vertexData[iOpt].floatArrays[iArray], swFloat,
								(maxCols+maxCells)*dims);
					vertexPtr[vIdx++] 
						= vertexData[iOpt].floatArrays[iArray];
					vMax++;
				}
				int tmpStart = maxCells;
				int icol;
				for(icol=0;icol<nSecs;icol++)
				{
					swInt colStart = colSAndC[2*icol];
					swInt colNum   = colSAndC[2*icol+1];
					DMA_IGet(&vertexData[iOpt].floatArrays
							[iArray][tmpStart*dims],
							&vertexData_m[iOpt].floatArrays
							[iArray][colStart*dims],
							colNum*dims*sizeof(swFloat), &status);
					dma_num++;
					tmpStart += colNum;
				}
				DMA_IGet(vertexData[iOpt].floatArrays[iArray],
							&vertexData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							nCells*dims*sizeof(swFloat), &status);
				dma_num++;
//if(_MYID==38) printf("%p,%f,%f\n",ldm_space_end,vertexData.floatArrays[iArray][maxCols+256],vertexData_m.floatArrays[iArray][30120]);
			}

			// selfConnData
			for(iArray=0;iArray<selfConnData[iOpt].fArrayNum;iArray++)
			{
				int dims = selfConnData[iOpt].fArrayDims[iArray];
				if(selfConnData[iOpt].fArrayInOut[iArray]==UPDATED)
				{
					selfConnData[iOpt].floatArrays[iArray] 
						= selfConnPtr[sIdx++];
					continue;
				} else if(sIdx<sMax)
				{
					selfConnData[iOpt].floatArrays[iArray] 
						= selfConnPtr[sIdx++];
				} else if(eIdx<eMax)
				{
					selfConnData[iOpt].floatArrays[iArray] 
						= edgePtr[eIdx++];
				} else
				{
					LDM_NEW(selfConnData[iOpt].floatArrays[iArray], swFloat,
								maxCells*dims);
					selfConnPtr[sIdx++] 
						= selfConnData[iOpt].floatArrays[iArray];
					sMax++;
				}
				DMA_IGet(selfConnData[iOpt].floatArrays[iArray],
							&selfConnData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							nCells*dims*sizeof(swFloat), &status);
				dma_num++;
			}

			// Get neighbor
			int iface,icell;
			for(iface=0;iface<nFaces;iface++)
			{
				owner[iface] -= cellStart;
				swInt col_global = neighbor[iface];
				swInt col_local  = 0;
				int isec1=0;
				while(col_global<colSAndC[2*isec1] ||
							col_global>=colSAndC[2*isec1+1]+colSAndC[2*isec1])
				{
					if(isec1<nSecs)
					{
						col_local += colSAndC[2*isec1+1];
						isec1++;
					} else
					{
						printf("Error at new neighbor!\n");
					}
				}
				neighbor[iface] = col_local+col_global-colSAndC[2*isec1]
					+ maxCells;
			}

			DMA_Wait(&status, dma_num);

			// Compute
			// Diagonal
//			frontEdgeData[iOpt].fArraySizes = 0;
//			backEdgeData[iOpt].fArraySizes = 0;
//			selfConnData[iOpt].fArraySizes  = nCells;
//			fun_slave(NULL,NULL,NULL,NULL,NULL,
//						owner,neighbor,&data_s[iOpt]);
//			swFloat *x = accessArray(&vertexData[iOpt], 0);
//			swFloat *b = accessArray(&vertexData[iOpt], 1);
//			swFloat *diag = accessArray(&selfConnData[iOpt], 0);
//				for(icell=0;icell<nCells;icell++)
//				{
////if(icell+cellStart==895) printf("%d,%d,%d,%f,%f,%f\n",_MYID,icell,cellStart,b[icell],diag[icell],x[icell]);
//					b[icell] = diag[icell]*x[icell];
//				}

			// FrontEdgeData
			frontEdgeData[iOpt].fArraySizes = nFaces;
			backEdgeData[iOpt].fArraySizes = nFaces;
			selfConnData[iOpt].fArraySizes  = nCells;
// if(_MYID==0)
// {
// if(_MYID==0) printf("%f\n",paraData.floatArrays[0][1]);
			fun_slave(NULL,NULL,NULL,NULL,&paraData,
						owner,neighbor,&data_s[iOpt]);
		// }

//			x = accessArray(&vertexData, 0);
//			b = accessArray(&vertexData, 1);
//			for(iface=0;iface<nFaces;iface++)
//			{
//				b[owner[iface]] += edge[iface]*x[neighbor[iface]];
//			}

			// vertexData(output)
			dma_num = 0;
			status = 0;
			for(iArray=0;iArray<vertexData[iOpt].fArrayNum;iArray++)
			{
				int dims = vertexData[iOpt].fArrayDims[iArray];
				if(vertexData[iOpt].fArrayInOut[iArray]==COPYIN) continue;
				DMA_IPut(&vertexData_m[iOpt].floatArrays
							[iArray][cellStart*dims],
							vertexData[iOpt].floatArrays[iArray],
							nCells*dims*sizeof(swFloat),&status);
				dma_num++;
			}
			DMA_Wait(&status, dma_num);
			int round;
			if(isec==0 && iOpt==0)
			{
			   	LDM_NEW(tmpData, swFloat, maxColNum*maxDims);
			}
			for(round=0;round<colRoundNum;round++)
			{
			//if(_MYID==0) printf("%f\n",vertexData_m[iOpt].floatArrays[0][260]);
				athread_syn(ARRAY_SCOPE,0xFFFF);
				if(round==colRound)
				{
					int icol,i_in_sec;
					int tmpStart = maxCells;
					dma_num = 0;
					status = 0;
					for(icol=0;icol<nSecs;icol++)
					{
						swInt colStart = colSAndC[2*icol];
						swInt colNum   = colSAndC[2*icol+1];
						for(iArray=0;iArray<vertexData[iOpt].fArrayNum;
									iArray++)
						{	
							int dims = vertexData[iOpt].fArrayDims[iArray];
							if(vertexData[iOpt].fArrayInOut[iArray]==COPYIN)
							{
								continue;
							}
							DMA_Get(tmpData,
								&vertexData_m[iOpt].floatArrays
								[iArray][colStart*dims],
								colNum*dims*sizeof(swFloat));
							swFloat *y = 
								&vertexData[iOpt].floatArrays
								[iArray][tmpStart*dims];
							for(i_in_sec=0;i_in_sec<colNum*dims;i_in_sec++)
							{
								y[i_in_sec] += tmpData[i_in_sec];
							}
							DMA_Put(&vertexData_m[iOpt].floatArrays
										[iArray][colStart*dims],
										y,
										colNum*dims*sizeof(swFloat));
//if(_MYID==0) printf("%f\n",vertexData_m[0].floatArrays[1][15]);
						}
						tmpStart += colNum;
					}
				}
			}
			athread_syn(ARRAY_SCOPE,0xFFFF);
		}
	}
}

void rowSubsectionIterator_array_slave(RSS_vertex2EdgePara *para)
{
// ************************************************************************
// public data structure
// ************************************************************************
#define MAX_LDM 62000
	RSS_vertex2EdgePara para_s = *para;
	const swInt faceNum        = para_s.edgeNumber;

	INIT_LDM_SPACE(MAX_LDM);

	// public data
	int iedge,dims,iDim;
	int iArray;
	int optNum = para_s.optNum;
	Arrays frontEdgeData_m[optNum], backEdgeData_m[optNum];
   	Arrays selfConnData_m[optNum], vertexData_m[optNum];
	Arrays frontEdgeData[optNum],   backEdgeData[optNum];
	Arrays selfConnData[optNum],   vertexData[optNum];
	Arrays paraData, paraData_m;
	FieldData data_s[optNum];
	int iOpt;
	for(iOpt=0;iOpt<optNum;iOpt++)
	{
		data_s[iOpt].frontEdgeData = &frontEdgeData[iOpt];
		data_s[iOpt].backEdgeData  = &backEdgeData[iOpt];
		data_s[iOpt].selfConnData  = &selfConnData[iOpt];
		data_s[iOpt].vertexData    = &vertexData[iOpt];
	}

	/// 获取参数及常量
	DMA_Get(&paraData_m,      para_s.paraData,      sizeof(Arrays));
	copyArrayWithoutData_new(&paraData_m, paraData,
				paraData_m.fArraySizes);
	for(iArray=0;iArray<paraData.fArrayNum;iArray++)
	{
		int length = paraData.fArraySizes; 
		LDM_NEW(paraData.floatArrays[iArray],swFloat,
				length);
		DMA_Get(paraData.floatArrays[iArray],
				paraData_m.floatArrays[iArray],
				length*sizeof(swFloat));
	}

	/// 暂时不提供融合算子功能
	iOpt = 0;

	FieldData data = *para_s.cOpt[iOpt].data;
	e2v_slaveFunPtr fun_slave = para_s.cOpt[iOpt].fun_slave;


	// 获取向量数组到从核，但不包括数组值
	DMA_Get(&frontEdgeData_m[iOpt], data.frontEdgeData, 
		sizeof(Arrays));
	copyArrayWithoutData_new(&frontEdgeData_m[iOpt], 
		frontEdgeData[iOpt], 0);

	int maxDims = 0;
	// FrontEdgeData
	for(iArray=0;iArray<frontEdgeData[iOpt].fArrayNum;iArray++)
	{
		maxDims += frontEdgeData[iOpt].fArrayDims[iArray];
	}
	swInt round = (faceNum*maxDims*sizeof(swFloat)-1)/(BLOCKNUM64K*MAX_LDM)+1;
	swInt faceNumInCPE = (0.96*MAX_LDM)/(sizeof(swFloat)*maxDims);
	swInt faceNumInRound = BLOCKNUM64K*faceNumInCPE;
	for(iArray=0;iArray<frontEdgeData[iOpt].fArrayNum;iArray++)
	{
		int dims = frontEdgeData[iOpt].fArrayDims[iArray];
		LDM_NEW(frontEdgeData[iOpt].floatArrays[iArray], 
			swFloat, faceNumInCPE*dims);
	}
// if(_MYID!=0) return;
// printf("%.16f,%p\n", frontEdgeData_m[0].floatArrays[0][76],data.frontEdgeData->floatArrays[0]);
	int iround,id,iface,j;
	for(iround=0;iround<round;iround++)
	{
		// for(id=0;id<64;id++)
		// {
		swInt faceStart = faceNumInRound*iround+_MYID*faceNumInCPE;
		// if(faceStart>=faceNum) break;
		faceNumInCPE = faceStart+faceNumInCPE > faceNum ? 
			faceNum-faceStart : faceNumInCPE;
// printf("%d,%d,%d,%d\n",iround,_MYID,faceStart,faceNum);
		for(iArray=0;iArray<frontEdgeData[iOpt].fArrayNum;iArray++)
		{
			int dims = frontEdgeData[iOpt].fArrayDims[iArray];
			if(frontEdgeData[iOpt].fArrayInOut[iArray]==COPYOUT) continue;
			DMA_Get(frontEdgeData[iOpt].floatArrays[iArray],
				&frontEdgeData_m[iOpt].floatArrays[iArray][faceStart*dims],
				faceNumInCPE*dims*sizeof(swFloat));
		}
		frontEdgeData[iOpt].fArraySizes = faceNumInCPE;
		// swInt *owner, *neighbor;
		fun_slave(NULL,NULL,NULL,NULL,&paraData,
					NULL,NULL,&data_s[iOpt]);
		// swFloat *gra = accessArray(&frontEdgeData[iOpt],0);
		// swInt dims_gra = getArrayDims(&frontEdgeData[iOpt],0);
		// swFloat *vol = accessArray(&frontEdgeData[iOpt],1);
		// for(iface=0;iface<faceNumInCPE;iface++)
		// {
		// 	for(j=0;j<18;j++)
		// 	{
		// 		gra[iface*dims_gra+j] /= vol[iface];
		// 	}
		// }
// if(iround==0 && id==0) printf("%f,%f\n", frontEdgeData_m[0].floatArrays[0][76], frontEdgeData[0].floatArrays[0][76]);
		for(iArray=0;iArray<frontEdgeData[iOpt].fArrayNum;iArray++)
		{
			int dims = frontEdgeData[iOpt].fArrayDims[iArray];
			if(frontEdgeData[iOpt].fArrayInOut[iArray]==COPYIN) continue;
			DMA_Put(&frontEdgeData_m[iOpt].floatArrays[iArray][faceStart*dims],
				frontEdgeData[iOpt].floatArrays[iArray],
				faceNumInCPE*dims*sizeof(swFloat));
		}
// if(iround==0 && id==0) printf("%f,%f\n", frontEdgeData_m[0].floatArrays[0][76], frontEdgeData[0].floatArrays[0][76]);
		// }
	}
#undef MAX_LDM
}
