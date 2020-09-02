#ifndef FILE_PRGRPCGC_SLAVE_C
#define FILE_PRGRPCGC_SLAVE_C
#include <slave.h>
#include <dma.h>
#include <unistd.h>
#include <simd.h>


#include "RlmpiSharedType.h"
#include "rlmpi.h"
#include "slaveUtils.h"


#if USE_DYNAMIC_MEM_INDICE==1
__thread_local_fix volatile __attribute__ ((aligned(32))) Pack *_sPacks;
__thread_local_fix volatile __attribute__ ((aligned(32))) Pack *_rPacks;
__thread_local_fix Pack *_sPacks_same_col;
__thread_local_fix Pack *_sPacks_same_row;

__thread_local volatile __attribute__ ((aligned(32))) int8LDM* _putr_schedules;
__thread_local volatile __attribute__ ((aligned(32))) int8LDM* _getrputc_schedules;
__thread_local volatile __attribute__ ((aligned(32))) int8LDM* _getc_schedules;


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _putr_schedules_same_row;
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _getr_schedules_same_row;


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _putc_schedules_same_col;
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _getc_schedules_same_col;

#else
__thread_local_fix Pack _sPacks[MaxNPackages];
__thread_local_fix Pack _rPacks[MaxNPackages];
__thread_local_fix Pack *_sPacks_same_col;
__thread_local_fix Pack *_sPacks_same_row;


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putr_schedules[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getrputc_schedules[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getc_schedules[MaxNCycle];

__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putr_schedules_same_row[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getr_schedules_same_row[MaxNCycle];


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putc_schedules_same_col[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getc_schedules_same_col[MaxNCycle];
#endif


__thread_local_fix volatile Table _table_ldm;
__thread_local_fix int _total_send_pcg;
__thread_local_fix int _total_recv_pcg;
__thread_local_fix int _nCycle;
__thread_local_fix int _nCycleSameCol;
__thread_local_fix volatile int _nCycleSameRow;

__thread_local int _get_reply, _put_reply;

inline void TransformPackage3(const Pack*sPacks, Pack*rPacks);
inline void TransformSameColumnPackage(const Pack*sPacks, Pack*rPacks);
inline void TransformSameRowPackage(const Pack*sPacks, Pack*rPacks);
inline void load_reg_mpi_init_data(Schedule * reg_data);
inline void REG_SIMD_GETR_PUTC();

inline void REG_SIMD_GETR_PUTC() {
    Pack var_;
    REG_SIMD_GETR(*(int256*) & var_);
    REG_SIMD_PUTC(*(int256*) & var_, ROW(var_.dst_id));
}

inline int largerest(int a, int b, int c) {
    int max = a;
    if (b > max) max = b;
    if (c > max) max = c;
    return max;
}

void sunway_check_memory_left(int *data) {

}
///////////////////Begin REG_MPI///////////////////////////
///////////////////////////////////////////////////////////////////////

inline void load_reg_mpi_init_data(Schedule * reg_data) {

    _total_send_pcg = _table_ldm.nPUTR + _table_ldm.nPutrSameRow + _table_ldm.nPutcSameCol;
    _total_recv_pcg = _table_ldm.nGETC + _table_ldm.nGetcSameCol + _table_ldm.nGetrSameRow;

#if USE_DYNAMIC_MEM==1
    int destory = reg_data->destroy;
    if (_total_send_pcg > 0 && !destory) {
        _sPacks = (Pack*) ldm_malloc(sizeof (Pack) * _total_send_pcg);
    } 
    
    if (_total_recv_pcg > 0 && !destory) {
        _rPacks = (Pack*) ldm_malloc(sizeof (Pack) * _total_recv_pcg);
    } 
#endif

    if (_total_send_pcg > 0) {
        _sPacks_same_row = &_sPacks[_table_ldm.nPUTR];
        _sPacks_same_col = &_sPacks[_table_ldm.nPUTR + _table_ldm.nPutrSameRow];
    }

    if (_total_send_pcg > 0) {
        _get_reply = 0;
        athread_get(PE_MODE, reg_data->package[_MYID],
                _sPacks, sizeof (Pack) * _total_send_pcg,
                &_get_reply, 0, 0, 0);
//        dma_wait(&_get_reply, 1);
		athread_wait(&_get_reply,1);
    }

}

void initTable(Schedule*reg_data) {

    _table_ldm.nGETC = reg_data->table[_MYID].nGETC;
    _table_ldm.nPUTR = reg_data->table[_MYID].nPUTR;
    _table_ldm.nGETR_PUTC = reg_data->table[_MYID].nGETR_PUTC;

    _nCycle = reg_data->nCycle;
    int destory = reg_data->destroy;

    int length = 32 * (sizeof (int8LDM) * _nCycle / 32 + 1);
    
    if (_nCycle > 0 && !destory) {
#if USE_DYNAMIC_MEM_INDICE==1
        _putr_schedules = (int8LDM*) ldm_malloc(length);
        _getrputc_schedules = (int8LDM*) ldm_malloc(length);
        _getc_schedules = (int8LDM*) ldm_malloc(length);
#endif
        _get_reply = 0;
        athread_get(PE_MODE, reg_data->putr_schedules[_MYID], _putr_schedules, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, reg_data->getrputc_schedules[_MYID], _getrputc_schedules, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, reg_data->getc_schedules[_MYID], _getc_schedules, length, &_get_reply, 0, 0, 0);
//        dma_wait(&_get_reply, 3);
		athread_wait(&_get_reply,3);

    } else if (_nCycle > 0 && destory) {
#if USE_DYNAMIC_MEM_INDICE==1
        ldm_free(_putr_schedules, length);
        ldm_free(_getrputc_schedules, length);
        ldm_free(_getc_schedules, length);
#endif
    }

    /////////////
    _nCycleSameRow = reg_data->nCycleSameRow;
    _table_ldm.nGetrSameRow = reg_data->table[_MYID].nGetrSameRow;
    _table_ldm.nPutrSameRow = reg_data->table[_MYID].nPutrSameRow;
    length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);
    if (_nCycleSameRow > 0 && !destory) {
#if USE_DYNAMIC_MEM_INDICE==1
        _putr_schedules_same_row = (int8LDM*) ldm_malloc(length);
        _getr_schedules_same_row = (int8LDM*) ldm_malloc(length);
#endif
        _get_reply = 0;
        athread_get(PE_MODE, reg_data->putr_schedules_same_row [_MYID], _putr_schedules_same_row, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, reg_data->getr_schedules_same_row[_MYID], _getr_schedules_same_row, length, &_get_reply, 0, 0, 0);
//        dma_wait(&_get_reply, 2);
		athread_wait(&_get_reply,2);

    } else if (_nCycleSameRow > 0 && destory) {
#if USE_DYNAMIC_MEM_INDICE==1
        ldm_free(_putr_schedules_same_row, length);
        ldm_free(_getr_schedules_same_row, length);
#endif
    }
    /////////////
    _nCycleSameCol = reg_data->nCycleSameCol;
    _table_ldm.nGetcSameCol = reg_data->table[_MYID].nGetcSameCol;
    _table_ldm.nPutcSameCol = reg_data->table[_MYID].nPutcSameCol;
    length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);
    if (_nCycleSameCol > 0 && !destory) {
#if USE_DYNAMIC_MEM_INDICE==1
        _putc_schedules_same_col = (int8LDM*) ldm_malloc(length);
        _getc_schedules_same_col = (int8LDM*) ldm_malloc(length);
#endif
        _get_reply = 0;
        athread_get(PE_MODE, reg_data->putc_schedules_same_col [_MYID], _putc_schedules_same_col, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, reg_data->getc_schedules_same_col[_MYID], _getc_schedules_same_col, length, &_get_reply, 0, 0, 0);
//        dma_wait(&_get_reply, 2);
		athread_wait(&_get_reply,2);

    } else if (_nCycleSameCol > 0 && destory) {
#if USE_DYNAMIC_MEM_INDICE==1
        ldm_free(_putc_schedules_same_col, length);
        ldm_free(_getc_schedules_same_col, length);
#endif
    }

    load_reg_mpi_init_data(reg_data);
    _total_send_pcg = _table_ldm.nPUTR + _table_ldm.nPutrSameRow + _table_ldm.nPutcSameCol;
    _total_recv_pcg = _table_ldm.nGETC + _table_ldm.nGetcSameCol + _table_ldm.nGetrSameRow;
	if(destory>0){
		if(_total_send_pcg>0){
			ldm_free(_sPacks,_total_send_pcg*sizeof(Pack));
		}
		if(_total_recv_pcg>0){
			ldm_free(_rPacks,_total_recv_pcg*sizeof(Pack));
		}
	}


    ALLSYN;
}

inline void TransformPackage3(const Pack*sPacks_, Pack*rPacks_) {

    int spr = 0, sgrpc = 0, sgc = 0;

    int cycle;
    for (cycle = 0; cycle < _nCycle; cycle++) {
        int nfpr;
        for (nfpr = 0; nfpr < _putr_schedules[cycle]; nfpr++) {
            REG_SIMD_PUTR(*(int256*) & sPacks_[spr], COL(sPacks_[spr].dst_id));
            spr++;
        }

        int nfgrpc;
        for (nfgrpc = 0; nfgrpc < _getrputc_schedules[cycle]; nfgrpc++) {
            REG_SIMD_GETR_PUTC();
        }

        int nfgc;
        for (nfgc = 0; nfgc < _getc_schedules[cycle]; nfgc++) {
            Pack tmp;
            REG_SIMD_GETC(*(int256*) & tmp);
            rPacks_[tmp.res_pos] = tmp;
        }
        ALLSYN;
    }

}

inline void TransformSameRowPackage(const Pack*sPacks, Pack*rPacks) {

    int spr = 0;

    int cycle;
    for (cycle = 0; cycle < _nCycleSameRow; cycle++) {

        int nfpr;
        for (nfpr = 0; nfpr < _putr_schedules_same_row[cycle]; nfpr++) {
            REG_SIMD_PUTR(*(int256*) & sPacks[spr], COL(sPacks[spr].dst_id));
            spr++;
        }

        int nfgr;
        for (nfgr = 0; nfgr < _getr_schedules_same_row[cycle]; nfgr++) {
            Pack tmp;
            REG_SIMD_GETR(*(int256*) & tmp);
            rPacks[tmp.res_pos] = tmp;
        }
    }
    ALLSYN;
}

inline void TransformSameColumnPackage(const Pack*sPacks, Pack*rPacks) {

    int spc = 0;

    int cycle;
    for (cycle = 0; cycle < _nCycleSameCol; cycle++) {

        int nfpc;
        for (nfpc = 0; nfpc < _putc_schedules_same_col[cycle]; nfpc++) {
            REG_SIMD_PUTC(*(int256*) & sPacks[spc], ROW(sPacks[spc].dst_id));
            spc++;
        }

        int nfgc;
        for (nfgc = 0; nfgc < _getc_schedules_same_col[cycle]; nfgc++) {
            Pack tmp;
            REG_SIMD_GETC(*(int256*) & tmp);
            rPacks[tmp.res_pos] = tmp;
        }
    }
    ALLSYN;
}

inline void sort_recv_package(Pack *rPacks_, int npack) {
    int j, i;
    for (i = 0; i < npack - 1; i++) {
        for (j = 0; j < npack - i - 1; j++) {
            if (rPacks_[j].src_id > rPacks_[j + 1].src_id) {
                Pack tmp;
                tmp = rPacks_[j];
                rPacks_[j] = rPacks_[j + 1];
                rPacks_[j + 1] = tmp;
            }
        }
    }
    ALLSYN;
}

void inline transform_data() {
    TransformPackage3(_sPacks, _rPacks);
    TransformSameRowPackage(_sPacks_same_row, _rPacks);
    TransformSameColumnPackage(_sPacks_same_col, _rPacks);
}
///////////////////////////////////////////////////////////////////
///////////////////End REG_MPI///////////////////////////
#endif
