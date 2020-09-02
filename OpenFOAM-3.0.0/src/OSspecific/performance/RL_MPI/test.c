/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <slave.h>
#include <dma.h>
#include <unistd.h>
#include <simd.h>
#include "rlmpi.h"

void test(Schedule* data) {
    //    total_send_pcg = table_ldm.nPUTR + table_ldm.nPutrSameRow + table_ldm.nPutcSameCol;
    //    int i = 0;
    //    for (i = 0; i < total_send_pcg; i++) {
    //        sPacks[i].data[0] = 0;
    //    }
    transform_data();
    int n;
    if (_MYID == 40) {
//        printf("%d\n", _total_recv_pcg);
        for (n = 0; n < _total_recv_pcg; n++) {
//            printf("%d\n", _rPacks[n].src_id);
        }
    }
}

void test_athread_get(double* in) {
    int get_reply = 0;
    double res[4];
    int i = 0;
    for (i = 0; i < 63 * 4 * 4; i++) {
        get_reply = 0;
        athread_get(PE_MODE, in, res, sizeof (double) *1, &get_reply, 0, 0, 0);
        dma_wait(&get_reply, 1);
    }
}
