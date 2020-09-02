#include "swMacro.h"
#include "slave.h"
#include <assert.h>
#include "rlmpi.h"
#include "slaveUtils.h"
#include "multiLevelBlockIterator.h"


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

#define initRegisterPacks(rlmpi_info) \
{  \
    _table_ldm.nGetcSkew = rlmpi_info->table[_MYID].nGetcSkew;\
    _table_ldm.nPutrSkew = rlmpi_info->table[_MYID].nPutrSkew;\
    _table_ldm.nGetrPutcSkew = rlmpi_info->table[_MYID].nGetrPutcSkew;\
    _nCycleSkew = rlmpi_info->nCycle;\
    int length = 32 * (sizeof (int8LDM) * _nCycleSkew / 32 + 1);\
    if (_nCycleSkew > 0) { \
		_putr_schedules_skew = (int8LDM*) ldm_malloc(length);\
		_getrputc_schedules_skew = (int8LDM*) ldm_malloc(length);\
		_getc_schedules_skew = (int8LDM*) ldm_malloc(length);\
		_get_reply = 0;\
		athread_get(PE_MODE, rlmpi_info->putr_schedules[_MYID], _putr_schedules_skew, length, &_get_reply, 0, 0, 0);\
		athread_get(PE_MODE, rlmpi_info->getrputc_schedules[_MYID], _getrputc_schedules_skew, length, &_get_reply, 0, 0, 0);\
		athread_get(PE_MODE, rlmpi_info->getc_schedules[_MYID], _getc_schedules_skew, length, &_get_reply, 0, 0, 0);\
		dma_wait(&_get_reply, 3);\
	} \
    _nCycleSameRow = rlmpi_info->nCycleSameRow;\
	_table_ldm.nGetrSameRow = rlmpi_info->table[_MYID].nGetrSameRow;\
	_table_ldm.nPutrSameRow = rlmpi_info->table[_MYID].nPutrSameRow;\
	length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);\
	if (_nCycleSameRow > 0) { \
		_putr_schedules_same_row = (int8LDM*) ldm_malloc(length);\
		_getr_schedules_same_row = (int8LDM*) ldm_malloc(length);\
		_get_reply = 0;\
		athread_get(PE_MODE, rlmpi_info->putr_schedules_same_row [_MYID], _putr_schedules_same_row, length, &_get_reply, 0, 0, 0); \
		athread_get(PE_MODE, rlmpi_info->getr_schedules_same_row[_MYID], _getr_schedules_same_row, length, &_get_reply, 0, 0, 0); \
		dma_wait(&_get_reply, 2);\
	} \
    _nCycleSameCol = rlmpi_info->nCycleSameCol;\
	_table_ldm.nGetcSameCol = rlmpi_info->table[_MYID].nGetcSameCol;\
	_table_ldm.nPutcSameCol = rlmpi_info->table[_MYID].nPutcSameCol;\
	length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);\
	if (_nCycleSameCol > 0) { \
		_putc_schedules_same_col = (int8LDM*) ldm_malloc(length);\
		_getc_schedules_same_col = (int8LDM*) ldm_malloc(length);\
		_get_reply = 0;\
		athread_get(PE_MODE, rlmpi_info->putc_schedules_same_col[_MYID], _putc_schedules_same_col, length, &_get_reply, 0, 0, 0);\
		athread_get(PE_MODE, rlmpi_info->getc_schedules_same_col[_MYID], _getc_schedules_same_col, length, &_get_reply, 0, 0, 0);\
		dma_wait(&_get_reply, 2);\
	} \
	load_rlmpi_data2(rlmpi_info);\
	ALLSYN;\
}

#define destroyRegisterPacks(rlmpi_info) \
{ \
	int length = 32 * (sizeof (int8LDM) * _nCycleSkew / 32 + 1);\
	if (_nCycleSkew > 0) {\
		ldm_free(_putr_schedules_skew, length);\
		ldm_free(_getrputc_schedules_skew, length);\
		ldm_free(_getc_schedules_skew, length);\
	}\
	_putr_schedules_skew = NULL;\
	_getrputc_schedules_skew = NULL;\
	_getc_schedules_skew = NULL;\
	\
    length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);\
	\
	if (length > 0) {\
		ldm_free(_putr_schedules_same_row, length);\
		ldm_free(_getr_schedules_same_row, length);\
	}\
	length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);\
	if (length > 0) {\
		ldm_free(_putc_schedules_same_col, length);\
		ldm_free(_getc_schedules_same_col, length);\
	}\
	if (_total_send_pcg > 0) {\
		ldm_free(_sPacks, _total_send_pcg * sizeof (Pack));\
	} \
	if (_total_recv_pcg > 0) {\
		ldm_free(_rPacks, _total_recv_pcg * sizeof (Pack));\
	}\
	ALLSYN;\
}

#define RLC_vertexData(vertexData_s) \
{ \
	swInt idx; \
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++) \
	{ \
		if(vertexData_s.fArrayInOut[iArray]==COPYOUT) continue; \
		if(vertexData_s.fArrayInOut[iArray]==UPDATED) continue; \
		dims = vertexData_s.fArrayDims[iArray]; \
		for(iDim=0;iDim<dims;iDim++) \
		{ \
			swFloat *x \
			    = vertexData_s.floatArrays[iArray]-cellStarts_s[0]*dims; \
			for(j=0;j<BLOCKNUM64K;j++) {sendIdx[j]=0;} \
			for(ipcg=0;ipcg<_total_send_pcg;ipcg++) \
			{ \
				if(_MYID>_sPacks[ipcg].dst_id) \
				{ \
					edgeIdx = recvList[recvIdx[_sPacks[ipcg].dst_id]] \
						+ sendIdx[_sPacks[ipcg].dst_id]; \
					for(j=0;j<_sPacks[ipcg].cva;j++) \
					{ \
						_sPacks[ipcg].data[j] \
							= x[neighbor[edgeIdx+j]*dims+iDim]; \
					} \
					sendIdx[_sPacks[ipcg].dst_id]+=_sPacks[ipcg].cva; \
				} \
			} \
			transform_data(); \
			recvNum = cellLen; \
			x += cellStarts_s[0]*dims; \
			for(ipcg=0;ipcg<_total_recv_pcg;ipcg++) \
			{ \
				if(_rPacks[ipcg].src_id>_MYID) \
				{ \
					for(j=0;j<_rPacks[ipcg].cva;j++) \
					{ \
						x[recvNum*dims+iDim] = _rPacks[ipcg].data[j]; \
						recvNum++; \
					} \
				} \
			} \
		} \
	} \
} 

#define RLC_vertexData_back(vertexData_s) \
{ \
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++) \
	{ \
		dims = vertexData_s.fArrayDims[iArray]; \
		if(vertexData_s.fArrayInOut[iArray]==COPYIN) continue; \
		if(vertexData_s.fArrayInOut[iArray]==UPDATED) continue; \
		swFloat* b = vertexData_s.floatArrays[iArray]; \
		for(iDim=0;iDim<dims;iDim++) \
		{ \
			for(j=0;j<BLOCKNUM64K;j++) \
			{ \
				sendIdx[j]=-blockStarts_s[6]+cellLen+cellStarts_s[0]; \
			} \
			for(ipcg=0;ipcg<_total_send_pcg;ipcg++) \
			{ \
				if(_MYID>=_sPacks[ipcg].dst_id) continue; \
				edgeIdx = blockStarts_s[4*(_sPacks[ipcg].dst_id-_MYID)+2] \
					+sendIdx[_sPacks[ipcg].dst_id]; \
				for(j=0;j<_sPacks[ipcg].cva;j++) \
				{ \
					_sPacks[ipcg].data[j] = b[(edgeIdx+j)*dims+iDim]; \
				} \
				sendIdx[_sPacks[ipcg].dst_id] += _sPacks[ipcg].cva; \
			} \
			transform_data(); \
			recvNum=0; \
			for(ipcg=0;ipcg<_total_recv_pcg;ipcg++) \
			{ \
				if(_MYID>_rPacks[ipcg].src_id) \
				{ \
					for(j=0;j<_rPacks[ipcg].cva;j++) \
					{ \
						b[neighbor[recvNum]*dims+iDim] \
							+= _rPacks[ipcg].data[j]; \
						recvNum++; \
					} \
				} else break; \
			} \
		} \
	} \
}

void multiLevelBlockIterator_e2v_slave(MLB_edge2VertexPara *para)
{
// ************************************************************************
// public data structure
// ************************************************************************
	MLB_edge2VertexPara para_s;
	DMA_Get(&para_s, para, sizeof(MLB_edge2VertexPara));

	int allocatable_ldm_size = 64000-(_total_send_pcg+_total_recv_pcg)*32;
//if(_MYID==32) printf("%d\n",allocatable_ldm_size);
	INIT_LDM_SPACE(54000);
	
	swInt *cellStarts_s, *blockStarts_s,
		  *owner_s, *neighbor_s;
	LDM_NEW(cellStarts_s,  swInt, 2);
	LDM_NEW(blockStarts_s, swInt, BLOCKNUM64K*4);

	swInt spIndex = para_s.spIndex;
	swInt mshBlockNum = para_s.mshBlockNum;
	swInt maxEdges = para_s.maxEdges;
	swInt maxCells = para_s.maxCells;
	swInt recvEdges = para_s.recvStarts[_MYID+spIndex*BLOCKNUM64K];
//	recvEdges = 0;

	swInt rowId = _MYID+spIndex*BLOCKNUM64K;
	swInt cpeBlockNum = mshBlockNum*BLOCKNUM64K;
	swInt startBlockIdx = rowId*(1+2*cpeBlockNum-rowId)/2;

	DMA_Get(blockStarts_s, &para_s.blockStarts[4*startBlockIdx],
				BLOCKNUM64K*4*sizeof(swInt));
	DMA_Get(cellStarts_s, &para_s.cellStarts[spIndex*BLOCKNUM64K+_MYID],
				2*sizeof(swInt));

	RlmpiInfo *rlmpi_info = &para_s.schedule_data[spIndex];
	initRegisterPacks(rlmpi_info);
#if 1

	int iArray;
	swInt *owner,*neighbor;
	int startIdx = blockStarts_s[2];
	int edgeLen = blockStarts_s[4*(BLOCKNUM64K-_MYID-1)+3]-startIdx;
	int localEdgeLen = blockStarts_s[3]-blockStarts_s[2];
	int cellLen = cellStarts_s[1]-cellStarts_s[0];
	int RL_length = cellLen+edgeLen-localEdgeLen;
	LDM_NEW(owner, swInt, edgeLen+recvEdges);
	LDM_NEW(neighbor, swInt, edgeLen+recvEdges);
	swInt *sendIdx,*recvIdx,*recvList,*sendIdxRev;
	LDM_NEW(sendIdx,  swInt, BLOCKNUM64K);
	LDM_NEW(recvIdx,  swInt, BLOCKNUM64K);
	LDM_NEW(recvList, swInt, BLOCKNUM64K);

	DMA_Get(&owner[recvEdges], &para_s.owner[startIdx],
				edgeLen*sizeof(swInt));
	DMA_Get(&neighbor[recvEdges], &para_s.neighbor[startIdx],
				edgeLen*sizeof(swInt));
	// Register communication: neighbor
	int ipcg,j,recvNum,edgeIdx;
	for(j=0;j<BLOCKNUM64K;j++) { sendIdx[j]=recvEdges-blockStarts_s[2]; }
	for(ipcg=0;ipcg<_total_send_pcg;ipcg++)
	{
		if(_MYID<_sPacks[ipcg].dst_id)
		{
			edgeIdx = blockStarts_s[4*(_sPacks[ipcg].dst_id-_MYID)+2]
				+ sendIdx[_sPacks[ipcg].dst_id];
			for(j=0;j<_sPacks[ipcg].cva;j++)
			{
				_sPacks[ipcg].data[j] = neighbor[edgeIdx+j];
			}
			sendIdx[_sPacks[ipcg].dst_id] += _sPacks[ipcg].cva;
		}
	}
	transform_data();
	recvNum=0;
	int k=0,lastId=0;
	for(j=0;j<BLOCKNUM64K;j++) { recvList[j]=0; }
	for(ipcg=0;ipcg<_total_recv_pcg;ipcg++)
	{
		if(_MYID>_rPacks[ipcg].src_id)
		{
			for(j=0;j<_rPacks[ipcg].cva;j++)
			{
				neighbor[recvNum++] = _rPacks[ipcg].data[j];
			}
			if(_rPacks[ipcg].src_id!=lastId)
			{
				k++;
				lastId=_rPacks[ipcg].src_id;
			}
			recvList[k+1]=recvNum;
			recvIdx[_rPacks[ipcg].src_id]=k;
		} else break;
	}

	// public data
	int iedge,dims,iDim;
	char* ldm_ptr = ldm_space_end;
	Arrays backEdgeData_c,frontEdgeData_c,vertexData_c,selfConnData_c;
	Arrays backEdgeData_m,frontEdgeData_m,vertexData_m,selfConnData_m;
	Arrays paraData, paraData_m;

	// localization
	swInt *p = &neighbor[recvEdges];
	for(iedge=localEdgeLen;iedge<edgeLen;iedge++)
	{
		p[iedge] = iedge-localEdgeLen+cellLen+cellStarts_s[0];
	}
	DMA_Get(&paraData_m,      para_s.paraData,      sizeof(Arrays));
	copyArrayWithoutData(&paraData_m, paraData,
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

// ************************************************************************
	int iOpt;
	int optNum = para_s.optNum;
	Arrays backEdgeData_s,frontEdgeData_s,vertexData_s,selfConnData_s;
	FieldData data_s;
	int eMax=0,vMax=0,sMax=0;
	swFloat* edgePtr[10];
	swFloat* vertexPtr[10];
	swFloat* selfConnPtr[10];

	char* ldm_tag = ldm_space_end;
	ldm_space_end += 1024;

	for(iOpt=0;iOpt<optNum;iOpt++) {

	ldm_ptr = ldm_tag;
	int eIdx=0,vIdx=0;

	FieldData *data = para_s.cOpt[iOpt].data_p;
	FieldData data_p;
	DMA_Get(&data_p, data, sizeof(FieldData));

//	DMA_Get(&backEdgeData_m,  data_p.backEdgeData,  sizeof(Arrays));
//	copyArrayWithoutData(&backEdgeData_m,  backEdgeData_c,
//				edgeLen);
//	for(iArray=0;iArray<backEdgeData_c.fArrayNum;iArray++)
//	{
//		if(backEdgeData_c.fArrayInOut[iArray]==UPDATED) continue;
//		int dims=backEdgeData_c.fArrayDims[iArray]; 
//		LDM_NEW(backEdgeData_c.floatArrays[iArray],swFloat,
//				edgeLen*dims);
//		DMA_Get(backEdgeData_c.floatArrays[iArray],
//				&backEdgeData_m.floatArrays[iArray][startIdx*dims],
//				edgeLen*dims*sizeof(swFloat));
//	}

	DMA_Get(&frontEdgeData_m, data_p.frontEdgeData, sizeof(Arrays));
	copyArrayWithoutData(&frontEdgeData_m, frontEdgeData_c,
				edgeLen);
	for(iArray=0;iArray<frontEdgeData_c.fArrayNum;iArray++)
	{
		int dims=frontEdgeData_c.fArrayDims[iArray]; 
		if(frontEdgeData_c.fArrayInOut[iArray]==UPDATED)
		{
			frontEdgeData_c.floatArrays[iArray] = edgePtr[eIdx++];
			continue;
		}
		if(eIdx<eMax)
		{
			frontEdgeData_c.floatArrays[iArray] = edgePtr[eIdx++];
		} else
		{
			LDM_NEW(frontEdgeData_c.floatArrays[iArray],swFloat,
					edgeLen*dims);
			edgePtr[eIdx++] = frontEdgeData_c.floatArrays[iArray];
			eMax++;
		}
		if(frontEdgeData_c.fArrayInOut[iArray]==COPYOUT) continue;
		DMA_Get(frontEdgeData_c.floatArrays[iArray],
				&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
				edgeLen*dims*sizeof(swFloat));
	}

//	DMA_Get(&selfConnData_m,  data_p.selfConnData,  sizeof(Arrays));
//	copyArrayWithoutData(&selfConnData_m, selfConnData_c, cellLen);
//	for(iArray=0;iArray<selfConnData_c.fArrayNum;iArray++)
//	{
//		int dims = getArrayDims(&selfConnData_c,iArray);
//		if(selfConnData_c.fArrayInOut[iArray]==UPDATED)
//		{
//			selfConnData_c.floatArrays[iArray] = selfConnPtr_p[iArray];
//			continue;
//		}
//		if(iArray<sMax_p)
//		{
//			selfConnData_c.floatArrays[iArray] = selfConnPtr_p[iArray];
//		} else
//		{
//			LDM_NEW(selfConnData_c.floatArrays[iArray],swFloat,
//					cellLen*dims);
//			selfConnPtr_p[iArray] = selfConnData_c.floatArrays[iArray];
//			sMax_p++;
//		}
//		DMA_Get(selfConnData_c.floatArrays[iArray],
//					&selfConnData_m.floatArrays
//					[iArray][cellStarts_s[0]*dims],
//					cellLen*dims*sizeof(swFloat));
//	}

	DMA_Get(&vertexData_m,    data_p.vertexData,    sizeof(Arrays));
	copyArrayWithoutData(&vertexData_m, vertexData_c, cellLen);
	for(iArray=0;iArray<vertexData_c.fArrayNum;iArray++)
	{
		int dims = getArrayDims(&vertexData_c,iArray);
		if(vertexData_c.fArrayInOut[iArray]==UPDATED)
		{
			vertexData_c.floatArrays[iArray] = vertexPtr[vIdx++];
			continue;
		}
		if(vIdx<vMax)
		{
			vertexData_c.floatArrays[iArray] = vertexPtr[vIdx++];
		} else if(eIdx < eMax)
		{
			vertexData_c.floatArrays[iArray] = edgePtr[eIdx++];
		} else
		{
			LDM_NEW(vertexData_c.floatArrays[iArray],swFloat,
						RL_length*dims);
			vertexPtr[vIdx++] = vertexData_c.floatArrays[iArray];
			vMax++;
		}
		DMA_Get(vertexData_c.floatArrays[iArray],
				&vertexData_m.floatArrays
				[iArray][cellStarts_s[0]*dims],
				cellLen*dims*sizeof(swFloat));
	}
	RLC_vertexData(vertexData_c);
	for(iArray=0;iArray<vertexData_c.fArrayNum;iArray++)
	{
		if(vertexData_c.fArrayInOut[iArray]==UPDATED) continue;
		dims = getArrayDims(&vertexData_c,iArray);
		if(vertexData_c.fArrayInOut[iArray]==COPYIN) continue;
		for(iedge=cellLen;iedge<RL_length;iedge++)
		{
			vertexData_c.floatArrays[iArray][iedge] = 0;
		}
	}	

	data = para_s.cOpt[iOpt].data;
	FieldData data_m;
	DMA_Get(&data_m, data, sizeof(FieldData));
	DMA_Get(&vertexData_m, data_m.vertexData,    sizeof(Arrays));
	DMA_Get(&selfConnData_m,data_m.selfConnData,  sizeof(Arrays));
	DMA_Get(&backEdgeData_m,data_m.backEdgeData,  sizeof(Arrays));
	DMA_Get(&frontEdgeData_m,data_m.frontEdgeData,sizeof(Arrays));
	copyArrayWithoutData(&vertexData_m, vertexData_s, cellLen);
	copyArrayWithoutData(&selfConnData_m, selfConnData_s, cellLen);
	copyArrayWithoutData(&backEdgeData_m,  backEdgeData_s, edgeLen);
	copyArrayWithoutData(&frontEdgeData_m, frontEdgeData_s, edgeLen);
	data_s.backEdgeData  = &backEdgeData_s;
	data_s.frontEdgeData = &frontEdgeData_s;
	data_s.selfConnData  = &selfConnData_s;
	data_s.vertexData    = &vertexData_s;

	e2v_slaveFunPtr fun_slave = para_s.cOpt[iOpt].fun_slave;
	// Construct array from master core without data copying
	backEdgeData_s.fArraySizes = edgeLen;
	frontEdgeData_s.fArraySizes = edgeLen;
	for(iArray=0;iArray<frontEdgeData_s.fArrayNum;iArray++)
	{
		int dims=frontEdgeData_s.fArrayDims[iArray]; 
		if(eIdx<eMax)
		{
			frontEdgeData_s.floatArrays[iArray] = edgePtr[eIdx];
		} else
		{
			LDM_NEW(frontEdgeData_s.floatArrays[iArray],swFloat,
					edgeLen*dims);
			edgePtr[eMax] = frontEdgeData_s.floatArrays[iArray];
			eMax++;
		}
		eIdx++;
		if(frontEdgeData_s.fArrayInOut[iArray]==COPYOUT) continue;
		DMA_Get(frontEdgeData_s.floatArrays[iArray],
				&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
				edgeLen*dims*sizeof(swFloat));
	}
	// Get vertex data through DMA
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vIdx<vMax)
		{
			vertexData_s.floatArrays[iArray] = vertexPtr[vIdx++];
		} else if(eIdx < eMax)
		{
			vertexData_s.floatArrays[iArray] = edgePtr[eIdx++];
		} else
		{
			LDM_NEW(vertexData_s.floatArrays[iArray],swFloat,
						RL_length*dims);
			vertexPtr[vIdx++] = vertexData_s.floatArrays[iArray];
			vMax++;
		}
		DMA_Get(vertexData_s.floatArrays[iArray],
			&vertexData_m.floatArrays[iArray][cellStarts_s[0]*dims],
			cellLen*dims*sizeof(swFloat));
	}
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vertexData_s.fArrayInOut[iArray]==COPYIN) continue;
		for(iedge=cellLen;iedge<RL_length;iedge++)
		{
			vertexData_s.floatArrays[iArray][iedge] = 0;
		}
	}

	// Get self connection data through DMA
	for(iArray=0;iArray<selfConnData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&selfConnData_s,iArray);
		if(iArray<sMax)
		{
			selfConnData_s.floatArrays[iArray] = selfConnPtr[iArray];
		} else
		{
			LDM_NEW(selfConnData_s.floatArrays[iArray],swFloat,
						cellLen*dims);
			selfConnPtr[iArray] = selfConnData_s.floatArrays[iArray];
			sMax++;
		}
		DMA_Get(selfConnData_s.floatArrays[iArray],
					&selfConnData_m.floatArrays
					[iArray][cellStarts_s[0]*dims],
					cellLen*dims*sizeof(swFloat));
	}
	// Register communication: vertex data(COPYIN)
	RLC_vertexData(vertexData_s);

	frontEdgeData_s.fArraySizes = 0;
	backEdgeData_s.fArraySizes = 0;
	selfConnData_s.fArraySizes = cellLen;
	frontEdgeData_c.fArraySizes = 0;
	backEdgeData_c.fArraySizes = 0;
	selfConnData_c.fArraySizes = cellLen;

	fun_slave(&backEdgeData_c, &frontEdgeData_c, &selfConnData_c,
				&vertexData_c, &paraData, 
				&owner[recvEdges], &neighbor[recvEdges], &data_s);

	// Compute Upper
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		vertexData_s.floatArrays[iArray] -= cellStarts_s[0]*dims;
	}
	for(iArray=0;iArray<vertexData_c.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_c,iArray);
		vertexData_c.floatArrays[iArray] -= cellStarts_s[0]*dims;
	}


	frontEdgeData_s.fArraySizes = edgeLen;
	backEdgeData_s.fArraySizes = 0;
	selfConnData_s.fArraySizes = 0;
	frontEdgeData_c.fArraySizes = edgeLen;
	backEdgeData_c.fArraySizes = 0;
	selfConnData_c.fArraySizes = 0;

	fun_slave(&backEdgeData_c, &frontEdgeData_c, &selfConnData_c,
				&vertexData_c, &paraData, 
				&owner[recvEdges], &neighbor[recvEdges], &data_s);
//if(_MYID==0) printf("%f,%f\n",data_s.vertexData->floatArrays[1][0],frontEdgeData_c.floatArrays[0][0]);
//	swFloat *face_lam = accessArray(&frontEdgeData_c,0);
//	swFloat *vis      = accessArray(data_s.vertexData,0);
//	swFloat *visFace  = accessArray(data_s.frontEdgeData,0);
//	swFloat *face_D   = accessArray(data_s.frontEdgeData,1);
//	swFloat *phi = accessArray(data_s.vertexData,0);
//	swFloat *S   = accessArray(data_s.vertexData,1);
//	swFloat *paras = accessArray(&paraData, 0);
//	swFloat faceVal;
////	swFloat *upper = accessArray(&frontEdgeData_s,0);
//	dims = getArrayDims(&frontEdgeData_s,0);
////if(_MYID==25) printf("%f,%f\n",data->frontEdgeData->floatArrays[0][194620],upper[1148]);
////	b -= cellStarts_s[_MYID];
//if(_MYID==0) printf("%f,%f,%f\n",paraData.floatArrays[0][0],paraData.floatArrays[0][1],paraData.floatArrays[0][2]);
//	for(iedge=0;iedge<edgeLen;iedge++)		
//	{
//		faceVal = face_lam[iedge]*vis[neighbor[iedge+recvEdges]]+(1-face_lam[iedge])*vis[owner[iedge+recvEdges]];
//		faceVal *= face_D[iedge];
//		faceVal = paras[0]/paras[1]+(faceVal-paras[0])/paras[2];
//if(iedge+blockStarts_s[2]==568066) printf("upper: %d,%d,%d,%f,%f,%f\n",_MYID,owner[iedge+recvEdges],neighbor[iedge+recvEdges],faceVal,paras[0],paras[1]);
//	}

//if(_MYID==25) printf("%f,%f\n",data->frontEdgeData->floatArrays[0][194620],upper[1148]);
	// Output frontEdgeData through DMA
	int flag = 0;
	for(iArray=0;iArray<frontEdgeData_s.fArrayNum;iArray++)
	{
		int dims=frontEdgeData_s.fArrayDims[iArray]; 
		if(frontEdgeData_s.fArrayInOut[iArray]==COPYIN) continue;
		if(frontEdgeData_s.fArrayInOut[iArray]==UPDATED) continue;
		flag = 1;
		DMA_Put(&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
			frontEdgeData_s.floatArrays[iArray],
			edgeLen*dims*sizeof(swFloat));
	}
	for(iArray=0;iArray<frontEdgeData_c.fArrayNum;iArray++)
	{
		dims = getArrayDims(&frontEdgeData_c,iArray);
		if(frontEdgeData_c.fArrayInOut[iArray]==COPYIN) continue;
		if(frontEdgeData_c.fArrayInOut[iArray]==UPDATED) continue;
		flag = 1;
		DMA_Put(&para_s.cOpt[iOpt].data_p->frontEdgeData->floatArrays
				[iArray][startIdx*dims],
				frontEdgeData_c.floatArrays[iArray],
				edgeLen*dims*sizeof(swFloat));
	}
//	if(flag)
//	{
//		// 恢复现场
//		for(iArray=0;iArray<vertexData_c.fArrayNum;iArray++)
//		{
//			dims = getArrayDims(&vertexData_c,iArray);
//			vertexData_c.floatArrays[iArray] += cellStarts_s[0]*dims;
//		}
//
//	   	continue;
//	}

// ***********************************************************************
// Backward edge
// ***********************************************************************
	for(iArray=0;iArray<backEdgeData_s.fArrayNum;iArray++)
	{
		int dims=backEdgeData_s.fArrayDims[iArray]; 
		backEdgeData_s.floatArrays[iArray]
			= frontEdgeData_s.floatArrays[iArray];
		DMA_Get(backEdgeData_s.floatArrays[iArray],
				&backEdgeData_m.floatArrays[iArray][startIdx*dims],
				edgeLen*dims*sizeof(swFloat));
	}

	// Compute Lower
	backEdgeData_s.fArraySizes = edgeLen;
	selfConnData_s.fArraySizes = 0;
	frontEdgeData_s.fArraySizes = 0;
	backEdgeData_c.fArraySizes = edgeLen;
	selfConnData_c.fArraySizes = 0;
	frontEdgeData_c.fArraySizes = 0;

	fun_slave(&backEdgeData_c, &frontEdgeData_c, &selfConnData_c,
				&vertexData_c, &paraData,
				&owner[recvEdges], &neighbor[recvEdges], &data_s);

	//TODO
	RLC_vertexData_back(vertexData_s);
	RLC_vertexData_back(vertexData_c);

	// 恢复现场
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		vertexData_s.floatArrays[iArray] += cellStarts_s[0]*dims;
	}
	for(iArray=0;iArray<vertexData_c.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_c,iArray);
		vertexData_c.floatArrays[iArray] += cellStarts_s[0]*dims;
	}

//	b = accessArray(&tmpVertexData,0);
//	x = accessArray(&tmpVertexData,1);
//	swFloat *lower = accessArray(&backEdgeData_s,0);
//	dims = getArrayDims(&backEdgeData_s,0);
////	b -= cellStarts_s[_MYID];
//	swInt length = localEdgeLen+recvEdges;
//	for(iedge=0;iedge<length;iedge++)				
//	{
//		for(iDim=0;iDim<dims;iDim++)
//		{
//if(neighbor[iedge*dims+iDim]==114912) printf("lower: %d,%d,%f,%f,%f\n",iedge,blockStarts_s[2],b[neighbor[iedge*dims+iDim]],lower[iedge*dims+iDim],x[owner[iedge]*dims+iDim]);
////			printf("%d,%d,%d,%d,%d,%d,%d\n",iedge,edgeLen,neighbor[iedge],cellStarts_s[iblk+_MYID+1],maxCells,ownCellLen,neiCellLen);
//			b[(neighbor[iedge])*dims+iDim]
//				+= lower[iedge*dims+iDim]*x[owner[iedge]*dims+iDim];
//		}
//	}

	// Output vertexData through DMA
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vertexData_s.fArrayInOut[iArray]==COPYIN) continue;
		if(vertexData_s.fArrayInOut[iArray]==UPDATED) continue;
		DMA_Put(&vertexData_m.floatArrays[iArray][cellStarts_s[0]*dims],
				vertexData_s.floatArrays[iArray],
				cellLen*dims*sizeof(swFloat));
	}
	for(iArray=0;iArray<vertexData_c.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_c,iArray);
		if(vertexData_c.fArrayInOut[iArray]==COPYIN) continue;
		if(vertexData_c.fArrayInOut[iArray]==UPDATED) continue;
		DMA_Put(&para_s.cOpt[iOpt].data_p->vertexData->floatArrays
				[iArray][cellStarts_s[0]*dims],
				vertexData_c.floatArrays[iArray],
				cellLen*dims*sizeof(swFloat));
	}
//if(_MYID==0) printf("%d,%d\n",iOpt, RLC_num);
#endif
	}
	destroyRegisterPacks(rlmpi_info);
}

void multiLevelBlockIterator_array_slave(MLB_arrayPara *para)
{
// ************************************************************************
// public data structure
// ************************************************************************
	MLB_arrayPara para_s;
	DMA_Get(&para_s, para, sizeof(MLB_arrayPara));

	INIT_LDM_SPACE(62000);
	
	swInt *blockStarts_s;
	LDM_NEW(blockStarts_s, swInt, 2);

	swInt mshBlockNum = para_s.mshBlockNum;

	Arrays backEdgeData_s,frontEdgeData_s,vertexData_s,selfConnData_s;
	FieldData data_s;
	swFloat* edgePtr[30];
	swFloat* vertexPtr[10];
	swFloat* selfConnPtr[10];

	// public data
	int iedge,dims,iDim;
	int iArray;
	char* ldm_ptr = ldm_space_end;
	Arrays backEdgeData_c,frontEdgeData_c,vertexData_c,selfConnData_c;
	Arrays backEdgeData_m,frontEdgeData_m,vertexData_m,selfConnData_m;
	Arrays paraData, paraData_m;
	DMA_Get(&paraData_m,      para_s.paraData,      sizeof(Arrays));
//	DMA_Get(&frontEdgeData_m, para_s.frontEdgeData, sizeof(Arrays));
//	copyArrayWithoutData(&frontEdgeData_m, frontEdgeData_c,
//				edgeLen);
	copyArrayWithoutData(&paraData_m, paraData,
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

	swInt spIndex, iOpt;
	char *ldm_tag = ldm_space_end;
	int optNum = para_s.optNum;
	for(spIndex=0;spIndex<mshBlockNum;spIndex++)
	{


//	REINIT_LDM_SPACE();
	swInt rowId = _MYID+spIndex*BLOCKNUM64K;
	swInt cpeBlockNum = mshBlockNum*BLOCKNUM64K;
	swInt startBlockIdx = rowId*(1+2*cpeBlockNum-rowId)/2;
	swInt endBlockIdx = (rowId+1)*(1+2*cpeBlockNum-rowId-1)/2-1;

	blockStarts_s[0] = para_s.blockStarts[4*startBlockIdx+2];
	blockStarts_s[1] = para_s.blockStarts[4*endBlockIdx+3];


	int eMax=0,vMax=0,sMax=0;
	ldm_space_end = ldm_tag;


//	for(iArray=0;iArray<frontEdgeData_c.fArrayNum;iArray++)
//	{
//		int dims=frontEdgeData_c.fArrayDims[iArray]; 
//		LDM_NEW(frontEdgeData_c.floatArrays[iArray],swFloat,
//				edgeLen*dims);
//		if(frontEdgeData_c.fArrayInOut[iArray]==COPYOUT) continue;
//		DMA_Get(frontEdgeData_c.floatArrays[iArray],
//				&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
//				edgeLen*dims*sizeof(swFloat));
//	}
// ************************************************************************
	for(iOpt=0;iOpt<optNum;iOpt++)
	{

	ldm_ptr = ldm_space_end; 

	int startIdx = blockStarts_s[0];
	int edgeLen = blockStarts_s[1]-startIdx;

	FieldData *data = para_s.cOpt[iOpt].data;
	FieldData data_m;
	DMA_Get(&data_m, data, sizeof(FieldData));
	DMA_Get(&frontEdgeData_m,data_m.frontEdgeData,sizeof(Arrays));
	copyArrayWithoutData(&frontEdgeData_m, frontEdgeData_s, edgeLen);
	data_s.frontEdgeData = &frontEdgeData_s;

	array_slaveFunPtr fun_slave = para_s.cOpt[iOpt].fun_slave;

	int arrayNum = frontEdgeData_s.fArrayNum;
	int sliceNum = edgeLen*arrayNum*sizeof(swFloat)/(60000*0.98)+1;
	int edgesInSlice = edgeLen/sliceNum+1;
	int iSlice;
	for(iSlice=0;iSlice<sliceNum;iSlice++)
	{
	startIdx = iSlice*edgesInSlice+blockStarts_s[0];
	edgeLen = startIdx+edgesInSlice > blockStarts_s[1] ?
		blockStarts_s[1]-startIdx : edgesInSlice;
//if(_MYID==0) printf("%d,%d,%d,%d,%d,%d\n",iSlice, startIdx, edgeLen, sliceNum, blockStarts_s[0], blockStarts_s[1]);
	int eIdx=0,vOutIdx=0;
	// Construct array from master core without data copying
	frontEdgeData_s.fArraySizes = edgeLen;
	for(iArray=0;iArray<frontEdgeData_s.fArrayNum;iArray++)
	{
		int dims=frontEdgeData_s.fArrayDims[iArray]; 
		if(iArray<eMax)
		{
			frontEdgeData_s.floatArrays[iArray] = edgePtr[iArray];
		} else
		{
			LDM_NEW(frontEdgeData_s.floatArrays[iArray],swFloat,
					edgeLen*dims);
			edgePtr[iArray] = frontEdgeData_s.floatArrays[iArray];
			eMax++;
		}
		eIdx++;
		if(frontEdgeData_s.fArrayInOut[iArray]==COPYOUT) continue;
		DMA_Get(frontEdgeData_s.floatArrays[iArray],
				&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
				edgeLen*dims*sizeof(swFloat));
	}
	frontEdgeData_s.fArraySizes = edgeLen;
	fun_slave(&backEdgeData_c, &frontEdgeData_c, &selfConnData_c,
				&vertexData_c, &paraData, &data_s);
//if(_MYID==0) printf("%f,%f\n",data_s.vertexData->floatArrays[1][0],frontEdgeData_c.floatArrays[0][0]);
//	swFloat *face_lam = accessArray(&frontEdgeData_c,0);
//	swFloat *vis      = accessArray(data_s.vertexData,0);
//	swFloat *visFace  = accessArray(data_s.frontEdgeData,0);
//	swFloat *face_D   = accessArray(data_s.frontEdgeData,1);
//	swFloat *phi = accessArray(data_s.vertexData,0);
//	swFloat *S   = accessArray(data_s.vertexData,1);
//	swFloat *paras = accessArray(&paraData, 0);
//	swFloat faceVal;
////	swFloat *upper = accessArray(&frontEdgeData_s,0);
//	dims = getArrayDims(&frontEdgeData_s,0);
////if(_MYID==25) printf("%f,%f\n",data->frontEdgeData->floatArrays[0][194620],upper[1148]);
////	b -= cellStarts_s[_MYID];
//if(_MYID==0) printf("%f,%f,%f\n",paraData.floatArrays[0][0],paraData.floatArrays[0][1],paraData.floatArrays[0][2]);
//	for(iedge=0;iedge<edgeLen;iedge++)		
//	{
//		faceVal = face_lam[iedge]*vis[neighbor[iedge+recvEdges]]+(1-face_lam[iedge])*vis[owner[iedge+recvEdges]];
//		faceVal *= face_D[iedge];
//		faceVal = paras[0]/paras[1]+(faceVal-paras[0])/paras[2];
//if(iedge+blockStarts_s[2]==568066) printf("upper: %d,%d,%d,%f,%f,%f\n",_MYID,owner[iedge+recvEdges],neighbor[iedge+recvEdges],faceVal,paras[0],paras[1]);
//	}

//if(_MYID==25) printf("%f,%f\n",data->frontEdgeData->floatArrays[0][194620],upper[1148]);
	// Output frontEdgeData through DMA
	for(iArray=0;iArray<frontEdgeData_s.fArrayNum;iArray++)
	{
		int dims=frontEdgeData_s.fArrayDims[iArray]; 
		if(frontEdgeData_s.fArrayInOut[iArray]==COPYIN) continue;
//if(_MYID==0) printf("%p,%p,%d\n",&frontEdgeData_m.floatArrays[iArray][startIdx*dims],frontEdgeData_s.floatArrays[iArray], edgeLen);
		DMA_Put(&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
					frontEdgeData_s.floatArrays[iArray],
					edgeLen*dims*sizeof(swFloat));
	}
	}
	}
//	for(iArray=0;iArray<frontEdgeData_c.fArrayNum;iArray++)
//	{
//		dims = getArrayDims(&frontEdgeData_c,iArray);
//		if(frontEdgeData_c.fArrayInOut[iArray]==COPYIN) continue;
//		DMA_Put(&para_s.frontEdgeData->floatArrays
//				[iArray][startIdx*dims],
//				frontEdgeData_c.floatArrays[iArray],
//				edgeLen*dims*sizeof(swFloat));
//	}
	}
}


void multiLevelBlockIterator_v2e_slave(MLB_vertex2EdgePara *para)
{
	MLB_vertex2EdgePara para_s;
	DMA_Get(&para_s, para, sizeof(MLB_vertex2EdgePara));

	INIT_LDM_SPACE(54000);

	swInt spIndex = para_s.spIndex;

	swInt *blockStarts, *cellStarts, *firstEdgeVertice, *vertexNeighbor;

	swInt startIdx, startBlockIdx, edgeLen, cellLen, localEdgeLen;
	swInt cellStart_DMA, fullCellLen;
	swInt i,j;
	swInt rowId = _MYID+spIndex*BLOCKNUM64K;
	swInt cpeBlockNum = para_s.mshBlockNum*BLOCKNUM64K;
	startBlockIdx = rowId*cpeBlockNum+spIndex*BLOCKNUM64K;
	
	DMA_Status status = 0;
	LDM_NEW(blockStarts, swInt, BLOCKNUM64K*4);
	DMA_IGet(blockStarts, &para_s.blockStarts[4*startBlockIdx],
				BLOCKNUM64K*4*sizeof(swInt),&status);
	LDM_NEW(cellStarts, swInt, 2);
	DMA_IGet(cellStarts, &para_s.cellStarts[rowId], 2*sizeof(swInt), &status);
	DMA_Wait(&status,2);

	startIdx = blockStarts[2];
	edgeLen = blockStarts[4*(BLOCKNUM64K-1)+3]-startIdx;
	cellLen = cellStarts[1]-cellStarts[0];
	cellStart_DMA = blockStarts[4*_MYID+2]-startIdx;
	localEdgeLen = blockStarts[4*_MYID+3]-blockStarts[4*_MYID+2];
	fullCellLen = edgeLen-localEdgeLen+cellLen;

	RlmpiInfo *rlmpi_info = &para_s.schedule_data[spIndex];
	initRegisterPacks(rlmpi_info);

	LDM_NEW(firstEdgeVertice, swInt, edgeLen);
	DMA_Get(firstEdgeVertice, &para_s.firstEdgeVertice[startIdx], 
				edgeLen*sizeof(swInt));
	LDM_NEW(vertexNeighbor, swInt, edgeLen);
	DMA_Get(vertexNeighbor, &para_s.vertexNeighbor[startIdx],
				edgeLen*sizeof(swInt));
	// Localization of vertexNeighbor
	for(i=0;i<cellStart_DMA;i++) { vertexNeighbor[i] = i+cellStarts[0]; }
	for(i=cellStart_DMA;i<cellStart_DMA+localEdgeLen;i++)
	{
		vertexNeighbor[i] += cellStart_DMA;
	}
	for(i=cellStart_DMA+localEdgeLen;i<edgeLen;i++)
	{
		vertexNeighbor[i] = i-localEdgeLen+cellLen+cellStarts[0];
	}


	// Get Field data: frontEdgeData, backEdgeData....
	Arrays frontEdgeData, backEdgeData, selfConnData, vertexData;
	FieldData data_s = {&backEdgeData, &frontEdgeData, &selfConnData,
		&vertexData};
	// Field data in MPE: store the pointer and value
	Arrays frontEdgeData_m, backEdgeData_m, selfConnData_m, vertexData_m;
	Arrays paraData, paraData_m;
	DMA_Get(&paraData_m, para_s.paraData, sizeof(Arrays));
	copyArrayWithoutData_new(&paraData_m, paraData, paraData_m.fArraySizes);
	int iArray;
	for(iArray=0;iArray<paraData.fArrayNum;iArray++)
	{
		int length = paraData.fArraySizes;
		LDM_NEW(paraData.floatArrays[iArray], swFloat, length);
		DMA_Get(paraData.floatArrays[iArray], paraData_m.floatArrays[iArray],
					length*sizeof(swFloat));
	}

	int iOpt;
	int optNum = para_s.optNum;
	int eMax=0,vMax=0,sMax=0;
	swFloat *edgePtr[10];
	swFloat *vertexPtr[10];
	swFloat *selfConnPtr[10];

	for(iOpt=0;iOpt<optNum;iOpt++)
	{

	int eIdx=0,vIdx=0,sIdx=0;
	FieldData *data = para_s.cOpt[iOpt].data;
	e2v_slaveFunPtr fun_slave = para_s.cOpt[iOpt].fun_slave;
	FieldData data_m;
	DMA_Get(&data_m, data, sizeof(FieldData));

	status = 0;
	DMA_IGet(&frontEdgeData_m, data_m.frontEdgeData, sizeof(Arrays),&status);
	DMA_IGet(&backEdgeData_m, data_m.backEdgeData, sizeof(Arrays),&status);
	DMA_IGet(&vertexData_m, data_m.vertexData, sizeof(Arrays),&status);
	DMA_IGet(&selfConnData_m, data_m.selfConnData, sizeof(Arrays),&status);
	DMA_Wait(&status,4);
	copyArrayWithoutData_new(&selfConnData_m, selfConnData, edgeLen);
	copyArrayWithoutData_new(&vertexData_m, vertexData, edgeLen);
	copyArrayWithoutData_new(&frontEdgeData_m, frontEdgeData, edgeLen);
	copyArrayWithoutData_new(&backEdgeData_m, backEdgeData, edgeLen);

	// VertexData
	for(iArray=0;iArray<vertexData.fArrayNum;iArray++)
	{
		int dims = vertexData.fArrayDims[iArray];
		if(vertexData.fArrayInOut[iArray]==UPDATED)
		{
			vertexData.floatArrays[iArray] = vertexPtr[vIdx++];
			continue;
		} else if(vIdx<vMax)
		{
			vertexData.floatArrays[iArray] = vertexPtr[vIdx++];
		} else if(eIdx<eMax)
		{
			vertexData.floatArrays[iArray] = edgePtr[eIdx++];
		} else
		{
			LDM_NEW(vertexData.floatArrays[iArray], swFloat, fullCellLen*dims);
			vertexPtr[vIdx++] = vertexData.floatArrays[iArray];
			vMax++;
		}
		DMA_Get(&vertexData.floatArrays[iArray][cellStart_DMA*dims], 
					&vertexData_m.floatArrays[iArray][cellStarts[0]*dims],
					cellLen*dims*sizeof(swFloat));
	}


	status = 0;
	int dma_num = 0;
	// FrontEdgeData
	for(iArray=0;iArray<frontEdgeData.fArrayNum;iArray++)
	{
		int dims = frontEdgeData.fArrayDims[iArray];
		if(frontEdgeData.fArrayInOut[iArray]==UPDATED)
		{
			frontEdgeData.floatArrays[iArray] = edgePtr[eIdx++];
			continue;
		} else if(eIdx<eMax)
		{
			frontEdgeData.floatArrays[iArray] = edgePtr[eIdx++];
		} else
		{
			LDM_NEW(frontEdgeData.floatArrays[iArray], swFloat, edgeLen*dims);
			edgePtr[eIdx++] = frontEdgeData.floatArrays[iArray];
			eMax++;
		}
		if(frontEdgeData.fArrayInOut[iArray]==COPYOUT) continue;
		DMA_IGet(frontEdgeData.floatArrays[iArray], 
					&frontEdgeData_m.floatArrays[iArray][startIdx*dims],
					edgeLen*dims*sizeof(swFloat),&status);
		dma_num++;
	}

	// BackEdgeData
	for(iArray=0;iArray<backEdgeData.fArrayNum;iArray++)
	{
		int dims = backEdgeData.fArrayDims[iArray];
		LDM_NEW(backEdgeData.floatArrays[iArray], swFloat, edgeLen*dims);
		DMA_IGet(backEdgeData.floatArrays[iArray], 
					&backEdgeData_m.floatArrays[iArray][startIdx*dims],
					edgeLen*dims*sizeof(swFloat),&status);
		dma_num++;
	}

	// SelfConnData
	for(iArray=0;iArray<selfConnData.fArrayNum;iArray++)
	{
		int dims = selfConnData.fArrayDims[iArray];
		if(selfConnData.fArrayInOut[iArray]==UPDATED)
		{
			selfConnData.floatArrays[iArray] = selfConnPtr[sIdx++];
			continue;
		} else if(sIdx<sMax)
		{
			selfConnData.floatArrays[iArray] = selfConnPtr[sIdx++];
		} else if(eIdx<eMax)
		{
			selfConnData.floatArrays[iArray] = edgePtr[eIdx++];
		} else
		{
			LDM_NEW(selfConnData.floatArrays[iArray], swFloat, cellLen*dims);
			selfConnPtr[sIdx++] = selfConnData.floatArrays[iArray];
			sMax++;
		}
		DMA_IGet(selfConnData.floatArrays[iArray], 
					&selfConnData_m.floatArrays[iArray][cellStarts[0]*dims],
					cellLen*dims*sizeof(swFloat),&status);
		dma_num++;
	}


	// Register communication
	swInt idx, recvUpper, recvLower, dims, iDim, edgeIdx, ipcg; 
	swInt *sendIdx;
	LDM_NEW(sendIdx, swInt, BLOCKNUM64K);
	for(iArray=0;iArray<vertexData.fArrayNum;iArray++)
	{
		if(vertexData.fArrayInOut[iArray]==COPYOUT) continue;
		if(vertexData.fArrayInOut[iArray]==UPDATED) continue;
		dims = vertexData.fArrayDims[iArray];
		for(iDim=0;iDim<dims;iDim++)
		{
			swFloat *x = vertexData.floatArrays[iArray]
				- (cellStarts[0] - cellStart_DMA)*dims;
			for(j=0;j<BLOCKNUM64K;j++) { sendIdx[j]=-blockStarts[2]; }
			for(ipcg=0;ipcg<_total_send_pcg;ipcg++)
			{
				edgeIdx = blockStarts[4*_sPacks[ipcg].dst_id+2]
					+ sendIdx[_sPacks[ipcg].dst_id];
				for(j=0;j<_sPacks[ipcg].cva;j++)
				{
					_sPacks[ipcg].data[j]
						= x[firstEdgeVertice[edgeIdx+j]*dims+iDim];
				}
				sendIdx[_sPacks[ipcg].dst_id]+=_sPacks[ipcg].cva;
			}
			transform_data();
			recvLower = 0;
			recvUpper = cellStart_DMA+cellLen;
			x = vertexData.floatArrays[iArray];
			for(ipcg=0;ipcg<_total_recv_pcg;ipcg++)
			{
				if(_rPacks[ipcg].src_id>_MYID)
				{
					for(j=0;j<_rPacks[ipcg].cva;j++)
					{
						x[recvUpper*dims+iDim] = _rPacks[ipcg].data[j];
						recvUpper++;
					}
				} else
				{
					for(j=0;j<_rPacks[ipcg].cva;j++)
					{
						x[recvLower*dims+iDim] = _rPacks[ipcg].data[j];
						recvLower++;
					}
				}
			}
		}
	}

	for(iArray=0;iArray<vertexData.fArrayNum;iArray++)
	{
		int dims = vertexData.fArrayDims[iArray];
		vertexData.floatArrays[iArray] += cellStart_DMA*dims;
	}
	DMA_Wait(&status, dma_num);

	// Compute
	selfConnData.fArraySizes = cellLen;
	frontEdgeData.fArraySizes = 0;
	backEdgeData.fArraySizes = 0;
	fun_slave(NULL,NULL,NULL,NULL,NULL,
				firstEdgeVertice, vertexNeighbor, &data_s);

//	swFloat *b = accessArray(&vertexData, 1);
//	swFloat *x = accessArray(&vertexData, 0);
//	swFloat *diag = accessArray(&selfConnData, 0);
//	for(i=0;i<cellLen;i++)
//	{
//		b[i] += diag[i]*x[i];
//	}
	for(iArray=0;iArray<vertexData.fArrayNum;iArray++)
	{
		int dims = vertexData.fArrayDims[iArray];
		if(vertexData.fArrayInOut[iArray]==COPYIN ||
					vertexData.fArrayInOut[iArray]==UPDATED)
		{
			vertexData.floatArrays[iArray] 
				-= (cellStart_DMA+cellStarts[0])*dims;
		} else if(vertexData.fArrayInOut[iArray]==COPYOUT)
		{
			vertexData.floatArrays[iArray] -= cellStarts[0]*dims;
		}
	}

	selfConnData.fArraySizes = 0;
	frontEdgeData.fArraySizes = edgeLen;
	backEdgeData.fArraySizes = 0;
	fun_slave(NULL,NULL,NULL,NULL,NULL,
				firstEdgeVertice, vertexNeighbor, &data_s);

//	swFloat *b = accessArray(&vertexData, 1);
//	swFloat *x = accessArray(&vertexData, 0);
//	swFloat *edgeData = accessArray(&frontEdgeData, 0);
//	swInt edgeNum = getArraySize(&frontEdgeData);
//	for(i=0;i<edgeLen;i++)
//	{
////if(firstEdgeVertice[i]==33971) printf("edge: %d,%f,%f,%f\n",i,b[firstEdgeVertice[i]],edgeData[i],x[vertexNeighbor[i]]);
//		b[firstEdgeVertice[i]] += edgeData[i]*x[vertexNeighbor[i]];
//	}

	// Output vertexData through DMA
	for(iArray=0;iArray<vertexData.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData,iArray);
		if(vertexData.fArrayInOut[iArray]==COPYIN) continue;
		if(vertexData.fArrayInOut[iArray]==UPDATED) continue;
		DMA_Put(&vertexData_m.floatArrays[iArray][cellStarts[0]*dims],
				&vertexData.floatArrays[iArray][cellStarts[0]*dims],
				cellLen*dims*sizeof(swFloat));
	}
	for(iArray=0;iArray<vertexData.fArrayNum;iArray++)
	{
		int dims = vertexData.fArrayDims[iArray];
		if(vertexData.fArrayInOut[iArray]==COPYIN ||
					vertexData.fArrayInOut[iArray]==UPDATED)
		{
			vertexData.floatArrays[iArray] += cellStarts[0]*dims;
		} else if(vertexData.fArrayInOut[iArray]==COPYOUT)
		{
			vertexData.floatArrays[iArray] 
				+= (cellStarts[0]-cellStart_DMA)*dims;
		}
	}

	}

	destroyRegisterPacks(rlmpi_info);

}
