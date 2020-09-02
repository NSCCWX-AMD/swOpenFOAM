#ifndef FILE_RLMPI_C
#define FILE_RLMPI_C
#include <slave.h>
#include <dma.h>
#include <unistd.h>
#include <simd.h>

#include "RlmpiShared.h"
#include "rlmpi.h"

#if USE_DYNAMIC_MEM_INDICE==1
__thread_local_fix volatile __attribute__ ((aligned(32))) Pack *_sPacks;
__thread_local_fix volatile __attribute__ ((aligned(32))) Pack *_rPacks;
__thread_local_fix Pack *_sPacks_same_col;
__thread_local_fix Pack *_sPacks_same_row;

__thread_local volatile __attribute__ ((aligned(32))) int8LDM* _putr_schedules_skew;
__thread_local volatile __attribute__ ((aligned(32))) int8LDM* _getrputc_schedules_skew;
__thread_local volatile __attribute__ ((aligned(32))) int8LDM* _getc_schedules_skew;


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _putr_schedules_same_row;
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _getr_schedules_same_row;


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _putc_schedules_same_col;
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM* _getc_schedules_same_col;

#else
__thread_local_fix Pack _sPacks[MaxNPackages];
__thread_local_fix Pack _rPacks[MaxNPackages];
__thread_local_fix Pack *_sPacks_same_col;
__thread_local_fix Pack *_sPacks_same_row;


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putr_schedules_skew[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getrputc_schedules_skew[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getc_schedules_skew[MaxNCycle];

__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putr_schedules_same_row[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getr_schedules_same_row[MaxNCycle];


__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _putc_schedules_same_col[MaxNCycle];
__thread_local_fix volatile __attribute__ ((aligned(32))) int8LDM _getc_schedules_same_col[MaxNCycle];
#endif


__thread_local_fix volatile Table _table_ldm;
__thread_local_fix int _total_send_pcg;
__thread_local_fix int _total_recv_pcg;
__thread_local_fix int _nCycleSkew;
__thread_local_fix int _nCycleSameCol;
__thread_local_fix volatile int _nCycleSameRow;

__thread_local int _get_reply;

inline void TransformPackageSkew(const Pack*sPacks, Pack*rPacks);
inline void TransformSameColumnPackage(const Pack*sPacks, Pack*rPacks);
inline void TransformSameRowPackage(const Pack*sPacks, Pack*rPacks);
inline void load_rlmpi_data(RlmpiInfo * rlmpi_info);
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

inline int Ceil(int a, int b) {
    return (a / b) + ((a % b) != 0);
}

inline int get_aligned_length(int len) {
    return Ceil(len, 32.0)*32;
}
///////////////////Begin RLMPI///////////////////////////
///////////////////////////////////////////////////////////////////////

inline void load_rlmpi_data(RlmpiInfo * rlmpi_info) {

    _total_send_pcg = _table_ldm.nPutrSkew + _table_ldm.nPutrSameRow + _table_ldm.nPutcSameCol;
    _total_recv_pcg = _table_ldm.nGetcSkew + _table_ldm.nGetcSameCol + _table_ldm.nGetrSameRow;

#if USE_DYNAMIC_MEM==1
    if (_total_send_pcg > 0) {
        _sPacks = (Pack*) ldm_malloc(sizeof (Pack) * _total_send_pcg);
    }

    if (_total_recv_pcg > 0) {
        _rPacks = (Pack*) ldm_malloc(sizeof (Pack) * _total_recv_pcg);
    }
#endif

    if (_total_send_pcg > 0) {
        _sPacks_same_row = &_sPacks[_table_ldm.nPutrSkew];
        _sPacks_same_col = &_sPacks[_table_ldm.nPutrSkew + _table_ldm.nPutrSameRow];
    }

    if (_total_send_pcg > 0) {
        _get_reply = 0;
        athread_get(PE_MODE, rlmpi_info->package[_MYID],
                _sPacks, sizeof (Pack) * _total_send_pcg,
                &_get_reply, 0, 0, 0);
        dma_wait(&_get_reply, 1);
    }
}

inline void load_rlmpi_data2(RlmpiInfo * rlmpi_info) {

    _total_send_pcg = _table_ldm.nPutrSkew + _table_ldm.nPutrSameRow + _table_ldm.nPutcSameCol;
    _total_recv_pcg = _table_ldm.nGetcSkew + _table_ldm.nGetcSameCol + _table_ldm.nGetrSameRow;

#if USE_DYNAMIC_MEM==1
    if (_total_send_pcg > 0) {
        _sPacks = (Pack*) ldm_malloc(sizeof (Pack) * _total_send_pcg);
    }

    if (_total_recv_pcg > 0) {
        _rPacks = (Pack*) ldm_malloc(sizeof (Pack) * _total_recv_pcg);
    }
#endif

    if (_total_send_pcg > 0) {
        _sPacks_same_row = &_sPacks[_table_ldm.nPutrSkew];
        _sPacks_same_col = &_sPacks[_table_ldm.nPutrSkew + _table_ldm.nPutrSameRow];
    }

    if (_total_send_pcg > 0) {

        int len1 = get_aligned_length(sizeof (int8LDM) * _total_send_pcg);
        int len2 = get_aligned_length(sizeof (int16LDM) * _total_send_pcg);
        
        int8LDM* srcId_list = (int8LDM*) ldm_malloc(len1);
        int8LDM* dstId_list = (int8LDM*) ldm_malloc(len1);

        int16LDM* resPos_list = (int16LDM*) ldm_malloc(len2);
        int16LDM* cvb_list = (int16LDM*) ldm_malloc(len2);
        int16LDM* cva_list = (int16LDM*) ldm_malloc(len2);
        _get_reply = 0;
        athread_get(PE_MODE, rlmpi_info->srcId_list[_MYID],
                srcId_list, len1,
                &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->dstId_list[_MYID],
                dstId_list, len1,
                &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->resPos_list[_MYID],
                resPos_list, len2,
                &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->cva_list[_MYID],
                cva_list, len2,
                &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->cvb_list[_MYID],
                cvb_list, len2,
                &_get_reply, 0, 0, 0);
        dma_wait(&_get_reply, 5);
        int n;
        for (n = 0; n < _total_send_pcg; n++) {
            _sPacks[n].dst_id = dstId_list[n];
            _sPacks[n].src_id = srcId_list[n];
            _sPacks[n].res_pos = resPos_list[n];
            _sPacks[n].cva = cva_list[n];
            _sPacks[n].cvb = cvb_list[n];
        }

        ldm_free(srcId_list, len1);
        ldm_free(dstId_list, len1);
        ldm_free(resPos_list, len2);
        ldm_free(cva_list, len2);
        ldm_free(cvb_list, len2);

    }
}

void destroyRlmpiInfo(RlmpiInfo*rlmpi_info) {

    int length = 32 * (sizeof (int8LDM) * _nCycleSkew / 32 + 1);

    if (_nCycleSkew > 0) {
#if USE_DYNAMIC_MEM_INDICE==1
        ldm_free(_putr_schedules_skew, length);
        ldm_free(_getrputc_schedules_skew, length);
        ldm_free(_getc_schedules_skew, length);
#endif
    }
    _putr_schedules_skew = NULL;
    _getrputc_schedules_skew = NULL;
    _getc_schedules_skew = NULL;
    /////////////

    length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);

#if USE_DYNAMIC_MEM_INDICE==1
    if (length > 0) {
        ldm_free(_putr_schedules_same_row, length);
        ldm_free(_getr_schedules_same_row, length);
    }
#endif

    /////////////
    length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);

    if (length > 0) {
#if USE_DYNAMIC_MEM_INDICE==1
        ldm_free(_putc_schedules_same_col, length);
        ldm_free(_getc_schedules_same_col, length);
#endif
    }
#if USE_DYNAMIC_MEM_INDICE==1
    if (_total_send_pcg > 0) {
        ldm_free(_sPacks, _total_send_pcg * sizeof (Pack));
    }
    if (_total_recv_pcg > 0) {
        ldm_free(_rPacks, _total_recv_pcg * sizeof (Pack));
    }
#endif
    ALLSYN;
}

void initRlmpiInfo(RlmpiInfo*rlmpi_info) {

    _table_ldm.nGetcSkew = rlmpi_info->table[_MYID].nGetcSkew;
    _table_ldm.nPutrSkew = rlmpi_info->table[_MYID].nPutrSkew;
    _table_ldm.nGetrPutcSkew = rlmpi_info->table[_MYID].nGetrPutcSkew;

    _nCycleSkew = rlmpi_info->nCycle;

    int length = 32 * (sizeof (int8LDM) * _nCycleSkew / 32 + 1);

    if (_nCycleSkew > 0) {
#if USE_DYNAMIC_MEM_INDICE==1
        _putr_schedules_skew = (int8LDM*) ldm_malloc(length);
        _getrputc_schedules_skew = (int8LDM*) ldm_malloc(length);
        _getc_schedules_skew = (int8LDM*) ldm_malloc(length);
#endif
        _get_reply = 0;
        athread_get(PE_MODE, rlmpi_info->putr_schedules[_MYID], _putr_schedules_skew, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->getrputc_schedules[_MYID], _getrputc_schedules_skew, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->getc_schedules[_MYID], _getc_schedules_skew, length, &_get_reply, 0, 0, 0);
        dma_wait(&_get_reply, 3);
    }

    /////////////
    _nCycleSameRow = rlmpi_info->nCycleSameRow;
    _table_ldm.nGetrSameRow = rlmpi_info->table[_MYID].nGetrSameRow;
    _table_ldm.nPutrSameRow = rlmpi_info->table[_MYID].nPutrSameRow;
    length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);
    if (_nCycleSameRow > 0) {
#if USE_DYNAMIC_MEM_INDICE==1
        _putr_schedules_same_row = (int8LDM*) ldm_malloc(length);
        _getr_schedules_same_row = (int8LDM*) ldm_malloc(length);
#endif
        _get_reply = 0;
        athread_get(PE_MODE, rlmpi_info->putr_schedules_same_row [_MYID], _putr_schedules_same_row, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->getr_schedules_same_row[_MYID], _getr_schedules_same_row, length, &_get_reply, 0, 0, 0);
        dma_wait(&_get_reply, 2);
    }
    /////////////
    _nCycleSameCol = rlmpi_info->nCycleSameCol;
    _table_ldm.nGetcSameCol = rlmpi_info->table[_MYID].nGetcSameCol;
    _table_ldm.nPutcSameCol = rlmpi_info->table[_MYID].nPutcSameCol;
    length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);
    if (_nCycleSameCol > 0) {
#if USE_DYNAMIC_MEM_INDICE==1
        _putc_schedules_same_col = (int8LDM*) ldm_malloc(length);
        _getc_schedules_same_col = (int8LDM*) ldm_malloc(length);
#endif
        _get_reply = 0;
        athread_get(PE_MODE, rlmpi_info->putc_schedules_same_col[_MYID], _putc_schedules_same_col, length, &_get_reply, 0, 0, 0);
        athread_get(PE_MODE, rlmpi_info->getc_schedules_same_col[_MYID], _getc_schedules_same_col, length, &_get_reply, 0, 0, 0);
        dma_wait(&_get_reply, 2);
    }

    load_rlmpi_data2(rlmpi_info);
    ALLSYN;
}

inline void TransformPackageSkew(const Pack*sPacks_, Pack*rPacks_) {

    int spr = 0;

    int cycle;
    for (cycle = 0; cycle < _nCycleSkew; cycle++) {
        int nfpr;
        for (nfpr = 0; nfpr < _putr_schedules_skew[cycle]; nfpr++) {
            REG_SIMD_PUTR(*(int256*) & sPacks_[spr], COL(sPacks_[spr].dst_id));
            spr++;
        }

        int nfgrpc;
        for (nfgrpc = 0; nfgrpc < _getrputc_schedules_skew[cycle]; nfgrpc++) {
            REG_SIMD_GETR_PUTC();
        }

        int nfgc;
        for (nfgc = 0; nfgc < _getc_schedules_skew[cycle]; nfgc++) {
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
    TransformPackageSkew(_sPacks, _rPacks);
    TransformSameRowPackage(_sPacks_same_row, _rPacks);
    TransformSameColumnPackage(_sPacks_same_col, _rPacks);
}
///////////////////////////////////////////////////////////////////
///////////////////End REG_MPI///////////////////////////
#endif
