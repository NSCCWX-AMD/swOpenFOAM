/* 
 * File:   RlmpiInitializer.cpp
 * Author: Geng Chen
 * 
 * Created on November 30, 2017, 3:49 PM
 */
#ifndef FILE_RLMPI_INITIALLIZER_CXX
#define FILE_RLMPI_INITIALLIZER_CXX
#include "RlmpiInitializer.hxx"
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <utility> 
#include<string>
#include<fstream>
#include <sstream>

using namespace std;

template < class T >
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    //    os << "[";
    for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end() - 1; ++it) {
        os << "  " << *it << ", ";
    }
    os << *v.rbegin();
    os << ";";
    return os;
}

void generate_register_transform_table_skew(int(*dst_list)[64],
        int(*sendN)[64], Table *table) {
    int i, j;
    for (i = 0; i < 64; i++) {
        table[i].nGetcSkew = 0;
        table[i].nPutrSkew = 0;
        table[i].nGetrPutcSkew = 0;
    }

    int srcId, dstId;
    for (srcId = 0; srcId < 64; srcId++) {
        for (dstId = 0; dstId < 64; dstId++) {
            int dst = dst_list[srcId][dstId];

            if (dst >= 0) {
                int ns = sendN[srcId][dstId]; //number of sent package
                table[srcId].nPutrSkew += ns;

            }
        }
    }

    for (dstId = 0; dstId < 64; dstId++) {
        //search each source threads
        for (srcId = 0; srcId < 64; srcId++) {
            int dst = dst_list[srcId][dstId];

            if (dst == dstId) {
                int ns = sendN[srcId][dst]; //number of sent package
                table[dstId].nGetcSkew += ns;
            }
        }
    }

    for (int transId = 0; transId < 64; transId++) {

        int trans_row = ROW(transId);
        int trans_col = COL(transId);

        int clientId;
        //search in the row which I belong to
        for (clientId = 8 * trans_row; clientId < 8 * trans_row + 8; clientId++) {
            for (dstId = 0; dstId < 64; dstId++) {
                int dst = dst_list[clientId][dstId];
                if (dst >= 0) {
                    int dst_col = COL(dst);
                    if (dst_col == trans_col) {
                        int ns = sendN[clientId][dstId];
                        table[transId].nGetrPutcSkew += ns;
                    }
                }
            }
        }
    }
}

void generate_register_transform_table_same_row(const int(*dst_list)[64],
        const int(*sendN)[64], Table *table) {

    for (int i = 0; i < 64; i++) {
        table[i].nGetrSameRow = 0;
        table[i].nPutrSameRow = 0;
    }
    int srcId, dstId;

    for (srcId = 0; srcId < 64; srcId++) {
        for (dstId = 0; dstId < 64; dstId++) {
            int dst = dst_list[srcId][dstId];

            if (dst >= 0) {
                int ns = sendN[srcId][dstId]; //number of sent package
                table[srcId].nPutrSameRow += ns;
            }
        }
    }

    for (dstId = 0; dstId < 64; dstId++) {
        //search each source threads
        for (srcId = 0; srcId < 64; srcId++) {
            int dst = dst_list[srcId][dstId];

            if (dst == dstId) {
                int ns = sendN[srcId][dst]; //number of sent package
                table[dstId].nGetrSameRow += ns;
            }
        }
    }
    ///for test////
    //    for (srcId = 0; srcId < 64; srcId++) {
    //        cout << "Table srcId" << " " << table[srcId].nPutrSameRow << " " << table[srcId].nGetrSameRow << endl;
    //    }
}

void generate_register_transform_table_same_col(int(*dst_list)[64],
        int(*sendN)[64], Table *table) {
    int i, j;
    for (i = 0; i < 64; i++) {
        table[i].nGetcSameCol = 0;
        table[i].nPutcSameCol = 0;
    }
    int srcId, dstId;

    for (srcId = 0; srcId < 64; srcId++) {
        for (dstId = 0; dstId < 64; dstId++) {
            int dst = dst_list[srcId][dstId];

            if (dst >= 0) {
                int ns = sendN[srcId][dstId]; //number of sent packages
                table[srcId].nPutcSameCol += ns;
            }
        }
    }

    for (dstId = 0; dstId < 64; dstId++) {
        //search each source threads
        for (srcId = 0; srcId < 64; srcId++) {
            int dst = dst_list[srcId][dstId];

            if (dst == dstId) {
                int ns = sendN[srcId][dst]; //number of sent package
                table[dstId].nGetcSameCol += ns;
            }
        }
    }
    ///for test////
    //    for (srcId = 0; srcId < 64; srcId++) {
    //        cout << "Table srcId" << "  " << table[srcId].nPutcSameCol << " " << table[srcId].nGetcSameCol << endl;
    //    }
}

//constructor

RlmpiInitializer::RlmpiInitializer() {
    //    nThread = 64;
    //    table.resize(nThread);
}
//generate the packages for sending in each thread 

void RlmpiInitializer::generate_data() {
    srand(time(NULL));
    for (int myId = 0; myId < nThread; myId++) {
        //        if (myId == 0 || myId == 1 || myId == 2 || myId == 3 || myId == 4 || myId == 5 || myId == 6 || myId == 7)// for test
        //        if (myId == 0 || myId == 1 || myId == 9)
        //        if (myId == 0)
        for (int n = 0; n < destination_pool[myId].size(); n++) {

            int nSend = rand() % maxNPack; //how many packages will be sent to the destination
            //            int nSend = 10; //how many packages will be sent to the destination 
            //                if (n != 0) {
            //                    nSend = 0;
            //                }
            for (int i = 0; i < nSend; i++) {
                Pack pack;
                pack.src_id = myId;
                pack.dst_id = destination_pool[myId][n];
                //                pack.d_col = COL(pack.dst_id);
                //                pack.d_row = ROW(pack.dst_id);
                //                pack.s_col = COL(pack.src_id);
                //                pack.s_row = ROW(pack.src_id);
                pack.data[0] = double(myId) + double(i) / 1e2;
                pack.data[1] = double(myId) + double(i) / 1e3;
                pack.data[2] = double(myId) + double(i) / 1e4;
                int dstId = destination_pool[myId][n];
                packages[myId][dstId].push_back(pack);
            }
        }
    }
    //    cout << "Packages generated" << endl;
}
// first step is generate data

void RlmpiInitializer::generate_data_skew(const vector<vector < int8LDM> > &not_col_row_dst) {

    for (int myId = 0; myId < nThread; myId++) {

        for (int s = 0; s < not_col_row_dst[myId].size(); s++) {
            Pack pack;
            pack.src_id = myId;
            int dstId = not_col_row_dst[myId][s];
            pack.dst_id = dstId;
            //            pack.d_col = COL(pack.dst_id);
            //            pack.d_row = ROW(pack.dst_id);
            //            pack.s_col = COL(pack.src_id);
            //            pack.s_row = ROW(pack.src_id);

            //            pack.data[0] = double(myId) + double(s) / 1e2;
            //            pack.data[1] = double(myId) + double(s) / 1e3;
            //            pack.data[2] = double(myId) + double(s) / 1e4;
            pack.data[0] = float(dstId);
            pack.data[1] = float(dstId);
            pack.data[2] = float(dstId);
            pack.data[3] = float(dstId);
            pack.data[4] = float(dstId);
            pack.data[5] = float(dstId);

            packages[myId][dstId].push_back(pack);
        }
    }
    //    cout << "Packages generated" << endl;
}

void RlmpiInitializer::generate_table_skew() {
    for (int myId = 0; myId < nThread; myId++) {
        //        table[myId].maxOP = 0;
        table[myId].nGetcSkew = 0;
        table[myId].nGetrPutcSkew = 0;
        table[myId].nPutrSkew = 0;
        table[myId].nGetcSameCol = 0;
        table[myId].nPutcSameCol = 0;
        table[myId].nGetrSameRow = 0;
        table[myId].nGetrSameRow = 0;
    }

    int dst_list[64][64];
    int sendN_list[64][64];
    for (int myId = 0; myId < 64; myId++) {
        for (int dstId = 0; dstId < 64; dstId++) {
            sendN_list[myId][dstId] = 0;
            dst_list[myId][dstId] = -1;
        }
    }
    Table table_[64];
    for (int myId = 0; myId < nThread; myId++) {

        for (map<int, vector<Pack> >::iterator it = packages[myId].begin();
                it != packages[myId].end(); it++) {

            const vector<Pack> &packs = it->second;
            int dstId = it->first;
            int nSend = it->second.size();

            dst_list[myId][dstId] = dstId;
            sendN_list[myId][dstId] = nSend;
        }
    }

    generate_register_transform_table_skew(dst_list, sendN_list, table);

    //        for (int i = 0; i < nThread; i++) {
    //            cout << "ThreadId: " << i << "     nPr: " << table[i].nPUTR << "      nGRPC: " << table[i].nGETR_PUTC << "       nGC: " << table[i].nGETC << endl;
    //        }

    //    cout << "Generate table done" << endl;
}

vector<int> RlmpiInitializer::get_destination_pool(int myId) {
    vector<int> res;
    int my_row = ROW(myId);
    int my_col = COL(myId);
    for (int dstId = 0; dstId < nThread; dstId++) {
        int dst_row = ROW(dstId);
        int dst_col = COL(dstId);

        if (dst_col != my_col && dst_row != my_row) {
            res.push_back(dstId);
        }
    }
    if (res.size() != 49) {
        // cout << "error" << endl;
        abort();
    } else {
        return res;
    }
}

void RlmpiInitializer::reorder_packages_skew() {

    for (int myId = 0; myId < nThread; myId++) {

        vector<Pack> all_packs; // All packages sent of myId

        vector<int> sThread(nThread, 0); // counter of packages of dstId

        for (int n = 0; n < table[myId].nPutrSkew; n++) {
            for (map<int, vector<Pack> >::iterator it = packages[myId].begin();
                    it != packages[myId].end(); it++) {

                int dstId = it->first;
                int nSend = it->second.size(); // n packages sent to dstId
                const vector<Pack> &packs = it->second; //packages to dstId

                if (sThread[dstId] < nSend) {
                    all_packs.push_back(packs[sThread[dstId]++]); //see the excel in folder
                }
            }
        }
        res_packages[myId] = all_packs;
    }

    //    cout << "package reorder done!" << endl;
}

void RlmpiInitializer::reorder_packages_skew2() {
    vector<int> nextCycleStartPos(nThread);

    vector<vector<Pack> > res_packages2;

    res_packages2.resize(64);

    for (int myId = 0; myId < nThread; myId++) {
        res_packages2[myId].resize(res_packages[myId].size());
    }

    for (int i = 0; i < nThread; i++) {
        nextCycleStartPos[i] = 0;
    }
    int stillHasPackage = 1;
    while (stillHasPackage) {

        int s = 0;
        for (int myId = 0; myId < nThread; myId++) {
            if (nextCycleStartPos[myId] < res_packages[myId].size()) {
                s++;
            }
        }
        //stop if all packages has been sent
        if (s == 0) {
            stillHasPackage = 0;
            break;
        }


        for (int myId = 0; myId < nThread; myId++) {
            map<int, int> counter;
            for (int n = 0; n < nThread; n++) {
                int ind = nextCycleStartPos[myId] + n;
                if (ind < res_packages[myId].size()) {
                    int dst_id = res_packages[myId][ind].dst_id;
                    if (dst_id >= 0) {
                        counter[dst_id] = 1;
                    }
                }
            }

            vector<int> distance;
            for (int n = 0; n < counter.size(); n++) {
                int ind = nextCycleStartPos[myId] + n;
                if (ind < res_packages[myId].size()) {
                    distance.push_back(abs(res_packages[myId][ind].dst_id - myId));
                }
            }
            //obtain the index permutation after the sorting 
            vector<pair<int, int> > V;
            for (int n = 0; n < distance.size(); n++) {
                pair<int, int> P = make_pair(distance[n], n);
                V.push_back(P);
            }
            sort(V.begin(), V.end());
            //            std::reverse(V.begin(), V.end());
            ///exchange package
            for (int n = 0; n < counter.size(); n++) {
                int ind = nextCycleStartPos[myId] + V[n].second;
                int ind2 = nextCycleStartPos[myId] + n;
                res_packages2[myId][ind2] = res_packages[myId][ind];

                //                cout << ind2 << " " << ind << " ss" << endl;
            }

            nextCycleStartPos[myId] += counter.size();
        }

    }
    for (int myId = 0; myId < nThread; myId++) {
        for (int n = 0; n < res_packages2[myId].size(); n++) {
            res_packages[myId][n] = res_packages2[myId][n];
        }
    }
    //    cout << "reorder package2 done" << endl;

}

template <typename T>
std::string NumberToString(T Number) {
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}

void RlmpiInitializer::write_packages() {
    for (int myId = 0; myId < nThread; myId++) {
        string filename = ".csv";
        filename = NumberToString(myId) + filename;
        //        cout << filename << endl;
        ofstream fout(filename.c_str(), std::ofstream::out);
        for (int n = 0; n < res_packages[myId].size(); n++) {
            fout << res_packages[myId][n].dst_id << endl;
        }
    }
}

void RlmpiInitializer::write_schedule() {
    string filename = "schedules.csv";
    ofstream fout(filename.c_str(), std::ofstream::out);
    for (int cycle = 0; cycle < putr_schedules[0].size(); cycle++) {
        for (int myId = 0; myId < nThread; myId++) {

            fout << "Cycle: " << cycle << " ThreadId   " << myId <<
                    "  " << putr_schedules[myId][cycle] <<
                    "  " << getrputc_schedules[myId][cycle] <<
                    "  " << getc_schedules[myId][cycle] << endl;
        }
    }
}

//generate the putr getrputc and getc schedule for each thread

//generate the putr getrputc and getc schedule for each thread

void RlmpiInitializer::generate_schedule_skew() {

    putr_schedules.resize(nThread);
    getrputc_schedules.resize(nThread);
    getc_schedules.resize(nThread);

    vector<int> nextCycleStartPos(nThread, 0);

    vector<int> packagesLength(nThread);
    for (int myId = 0; myId < nThread; myId++) {
        packagesLength[myId] = res_packages[myId].size();
    }

    int cycle = 0;
    /////////////////////////////
    int stillHasPackage = 1;
    while (stillHasPackage) {

        int s = 0;
        for (int myId = 0; myId < nThread; myId++) {
            if (nextCycleStartPos[myId] < packagesLength[myId]) {
                s++;
            }
        }
        //stop if all packages has been sent
        if (s == 0) {
            stillHasPackage = 0;
            break;
        }

        //generate putr schedule
        for (int myId = 0; myId < nThread; myId++) {
            map<int, vector<int> > counter;
            int nFirstPr;
            //we can always do the first different 7 putr 
            for (int n = 0; n < 7; n++) {
                int ind = nextCycleStartPos[myId] + n;
                if (ind < res_packages[myId].size()) {
                    //                    int d_col = res_packages[myId][ind].d_col;
                    int d_col = COL(res_packages[myId][ind].dst_id);
                    int s = 1;
                    counter[d_col].push_back(s);
                }
            }
            nFirstPr = counter.size();

            putr_schedules[myId].push_back(nFirstPr);
        }


        ///////generate getrputc_schedules
        vector<vector<int> > getrputc_list;
        getrputc_list.resize(nThread);

        for (int myId = 0; myId < nThread; myId++) {

            //count how many getrputc will do in my thread 
            int myRow = ROW(myId);
            int myCol = COL(myId);
            getrputc_schedules[myId].push_back(0);
            //loop the threads in my row, find who will putr to myId
            for (int mySameRowId = myRow * 8; mySameRowId < myRow * 8 + 8; mySameRowId++) {
                if (mySameRowId != myId) {

                    int nFirstPr = putr_schedules[mySameRowId][cycle];

                    for (int n = 0; n < nFirstPr; n++) {
                        int ind = nextCycleStartPos[mySameRowId] + n;
                        //                        if (res_packages[mySameRowId][ind].d_col == myCol)
                        if (COL(res_packages[mySameRowId][ind].dst_id) == myCol) {
                            getrputc_schedules[myId][cycle]++;
                            //                            getrputc_list[myId].push_back(res_packages[mySameRowId][ind].d_row); //record which row myId will putc
                            getrputc_list[myId].push_back(ROW(res_packages[mySameRowId][ind].dst_id)); //record which row myId will putc
                        }
                    }
                }
            }
        }

        ///////generate getc schedules
        //TODO
        for (int myId = 0; myId < nThread; myId++) {

            //count how many getc will do in my thread 

            int myRow = ROW(myId);
            int myCol = COL(myId);
            getc_schedules[myId].push_back(0);
            for (int row = 0; row < 8; row++) {
                int mySameColId = myCol + row * 8;
                if (mySameColId != myId) {//loop the threads in my column 

                    int nFirstGrPc = getrputc_schedules[mySameColId][cycle]; //how many package this thread will transform

                    for (int n = 0; n < nFirstGrPc; n++) {
                        int d_row = getrputc_list[mySameColId][n];
                        // if the destination row is me
                        if (d_row == myRow) {
                            getc_schedules[myId][cycle]++;

                            //                            DISP2(cycle, getc_schedules[myId].size());
                            //                            printf("%d\n", getc_schedules[myId][cycle]);
                        }
                    }

                }
            }
        }

        for (int myId = 0; myId < nThread; myId++) {
            int nFirstPr = putr_schedules[myId][cycle];
            nextCycleStartPos[myId] += nFirstPr;
        }

        cycle++;
    }

    //    cout << "Generate schedules done!" << endl;
    //FOR TEST///////////
    for (int myId = 0; myId < nThread; myId++) {

        if (putr_schedules[myId].size() != getc_schedules[myId].size() ||
                putr_schedules[myId].size() != getrputc_schedules[myId].size() ||
                getc_schedules[myId].size() != getrputc_schedules[myId].size()) {
            cout << "error" << endl;
            abort();
        }
    }

    /////////////////////
    //    for (int i = 0; i < putr_schedules[0].size(); i++) {
    //        for (int myId = 0; myId < nThread; myId++) {
    //            cout << "Cycle: " << i << " ThreadId   " << myId << "  " << putr_schedules[myId][i] << "  " << getrputc_schedules[myId][i] << "  " << getc_schedules[myId][i] << endl;
    //        }
    //        cout << "\n\n\n" << endl;
    //    }
    /////////////////////
}

void RlmpiInitializer::generate_schedule_same_row() {



    vector<int> nextCycleStartPos(nThread);
    for (int myId = 0; myId < nThread; myId++) {
        nextCycleStartPos[myId] = 0;
    }
    vector<int> packagesLength(nThread);
    for (int myId = 0; myId < nThread; myId++) {
        packagesLength[myId] = res_packages_same_row[myId].size();
    }

    int cycle = 0;
    /////////////////////////////
    int stillHasPackage = 1;
    while (stillHasPackage) {

        int sq = 0;
        for (int myId = 0; myId < nThread; myId++) {
            if (nextCycleStartPos[myId] < packagesLength[myId]) {
                sq++;
            }
        }
        //stop if all packages has been sent
        if (sq == 0) {
            stillHasPackage = 0;
            break;
        }

        //generate putr schedule
        for (int myId = 0; myId < nThread; myId++) {
            map<int, vector<int> > counter;
            int nFirstPr;
            //we can always do the first different 7 putr 
            for (int n = 0; n < 7; n++) {
                int ind = nextCycleStartPos[myId] + n;
                if (ind < res_packages_same_row[myId].size()) {
                    int d_col = COL(res_packages_same_row[myId][ind].dst_id);
                    int s = 1;
                    counter[d_col].push_back(s);
                }
            }
            nFirstPr = counter.size();

            putr_schedules_same_row[myId].push_back(nFirstPr);
        }


        ///////generate getr_schedules
        vector<vector<int> > getrputc_list;
        getrputc_list.resize(nThread);

        for (int myId = 0; myId < nThread; myId++) {

            //count how many getrputc will do in my thread 
            int myRow = ROW(myId);
            int myCol = COL(myId);
            getr_schedules_same_row[myId].push_back(0);
            //loop the threads in my row, find who will putr to myId
            for (int mySameRowId = myRow * 8; mySameRowId < myRow * 8 + 8; mySameRowId++) {
                if (mySameRowId != myId) {

                    int nFirstPr = putr_schedules_same_row[mySameRowId][cycle];

                    for (int n = 0; n < nFirstPr; n++) {
                        int ind = nextCycleStartPos[mySameRowId] + n;
                        //                        if (res_packages_same_row[mySameRowId][ind].d_col == myCol)
                        if (COL(res_packages_same_row[mySameRowId][ind].dst_id) == myCol) {
                            getr_schedules_same_row[myId][cycle]++;
                            //                            getrputc_list[myId].push_back(res_packages_same_row[mySameRowId][ind].d_row); //record which row myId will putc
                            getrputc_list[myId].push_back(ROW(res_packages_same_row[mySameRowId][ind].dst_id)); //record which row myId will putc
                        }
                    }
                }
            }
        }





        for (int myId = 0; myId < nThread; myId++) {
            int nFirstPr = putr_schedules_same_row[myId][cycle];
            nextCycleStartPos[myId] += nFirstPr;
        }

        cycle++;
    }

    //    cout << "\nGenerate row schedules done!" << endl;
    //FOR TEST///////////
    for (int myId = 0; myId < nThread; myId++) {

        if (putr_schedules_same_row[myId].size() != getr_schedules_same_row[myId].size()) {
            cout << "error" << endl;
            abort();
        }
    }
    //    for (int cycle = 0; cycle < getr_schedules_same_row[0].size(); cycle++) {
    //        for (int myId = 0; myId < nThread; myId++) {
    //            cout << cycle << "  myId " << myId << " "
    //                    << putr_schedules_same_row[myId][cycle] << " "
    //                    << getr_schedules_same_row[myId][cycle] << endl;
    //        }
    //    }



}

void RlmpiInitializer::generate_data_same_row() {
    srand(time(NULL));
    for (int myId = 0; myId < nThread; myId++) {

        int myCol = COL(myId);
        int myRow = ROW(myId);
        for (int dstId = myRow * 8; dstId < myRow * 8 + 8; dstId++) {

            if (dstId != myId) {

                int nSend = rand() % maxNPack;
                //                int nSend = 100;

                for (int i = 0; i < nSend; i++) {
                    Pack pack;
                    pack.src_id = myId;
                    pack.dst_id = dstId;
                    //                    pack.d_col = COL(pack.dst_id);
                    //                    pack.d_row = ROW(pack.dst_id);
                    //                    pack.s_col = COL(pack.src_id);
                    //                    pack.s_row = ROW(pack.src_id);
                    //                    pack.data[0] = double(myId) + double(i) / 1e2;
                    //                    pack.data[1] = double(myId) + double(i) / 1e3;
                    //                    pack.data[2] = double(myId) + double(i) / 1e4;
                    packages[myId][dstId].push_back(pack);

                }
            }
        }

    }
    cout << "Packages generated" << endl;
}

void RlmpiInitializer::generate_data_same_col(const vector<vector <int8LDM> > &src2dst, const vector<vector <int8LDM> > &nData) {
    for (int myId = 0; myId < nThread; myId++) {
        same_col_packages[myId].clear();
        for (int s = 0; s < src2dst[myId].size(); s++) {
            Pack pack;
            pack.src_id = myId;
            int dstId = src2dst[myId][s];
            pack.dst_id = dstId;
            int d_col = COL(pack.dst_id);
            int d_row = ROW(pack.dst_id);
            int s_col = COL(pack.src_id);
            int s_row = ROW(pack.src_id);
            pack.cva = nData[myId][s];
            if (d_row == s_row || s_col != d_col || pack.dst_id == pack.src_id) {
                DISP("error");
                abort();
            }

            same_col_packages[myId][dstId].push_back(pack);
        }
    }
    //    cout << "same column packages generated" << endl;
}

void RlmpiInitializer::generate_data_same_row(const vector<vector <int8LDM> > &src2dst, const vector<vector <int8LDM> > &nData) {
    for (int myId = 0; myId < nThread; myId++) {
        same_row_packages[myId].clear();
        for (int s = 0; s < src2dst[myId].size(); s++) {
            Pack pack;
            pack.src_id = myId;
            int dstId = src2dst[myId][s];
            pack.dst_id = dstId;
            int d_col = COL(pack.dst_id);
            int d_row = ROW(pack.dst_id);
            int s_col = COL(pack.src_id);
            int s_row = ROW(pack.src_id);
            pack.cva = nData[myId][s];

            if (d_row != s_row || s_col == d_col || pack.dst_id == pack.src_id) {
                DISP("error");
                abort();
            }

            same_row_packages[myId][dstId].push_back(pack);
        }
    }
    //    cout << "same row packages generated" << endl;
}

void RlmpiInitializer::reorder_packages_same_row() {


    for (int myId = 0; myId < nThread; myId++) {

        vector<Pack> all_packs; // All packages sent of myId

        vector<int> sThread(nThread, 0); // counter of packages of dstId

        for (int n = 0; n < table[myId].nPutrSameRow; n++) {
            for (map<int, vector<Pack> >::iterator it = same_row_packages[myId].begin();
                    it != same_row_packages[myId].end(); it++) {

                int dstId = it->first;
                int nSend = it->second.size(); // n packages sent to dstId
                const vector<Pack> &packs = it->second; //packages to dstId

                if (sThread[dstId] < nSend) {
                    all_packs.push_back(packs[sThread[dstId]++]); //see the excel in folder
                }
            }
        }
        res_packages_same_row[myId] = all_packs;
    }

    //    cout << "same row package reorder done!" << endl;
}

void RlmpiInitializer::generate_table_same_row() {
    for (int myId = 0; myId < nThread; myId++) {
        table[myId].nGetrSameRow = 0;
        table[myId].nGetrSameRow = 0;
    }

    int dst_list[64][64];
    int sendN_list[64][64];
    for (int myId = 0; myId < 64; myId++) {
        for (int dstId = 0; dstId < 64; dstId++) {
            sendN_list[myId][dstId] = 0;
            dst_list[myId][dstId] = -1;
        }
    }

    for (int myId = 0; myId < nThread; myId++) {
        for (map<int, vector<Pack> >::iterator it = same_row_packages[myId].begin();
                it != same_row_packages[myId].end(); it++) {
            int dstId = it->first;
            int nSend = it->second.size();
            dst_list[myId][dstId] = dstId;
            sendN_list[myId][dstId] = nSend;
        }
    }
    generate_register_transform_table_same_row(dst_list, sendN_list, table);
    //    cout << "Generate same row table done" << endl;
}

void RlmpiInitializer::generate_schedule_same_col() {

    vector<int> nextCycleStartPos(nThread, 0);

    vector<int> packagesLength(nThread);
    for (int myId = 0; myId < nThread; myId++) {
        packagesLength[myId] = res_packages_same_col[myId].size();
    }

    int cycle = 0;
    /////////////////////////////
    int stillHasPackage = 1;
    while (stillHasPackage) {

        int sq = 0;
        for (int myId = 0; myId < nThread; myId++) {
            if (nextCycleStartPos[myId] < packagesLength[myId]) {
                sq++;
            }
        }
        //stop if all packages has been sent
        if (sq == 0) {
            stillHasPackage = 0;
            break;
        }

        //generate putc schedule
        for (int myId = 0; myId < nThread; myId++) {
            map<int, vector<int> > counter;
            int nFirstPr;
            //we can always do the first different 7 putc 
            for (int n = 0; n < 7; n++) {
                int ind = nextCycleStartPos[myId] + n;
                if (ind < res_packages_same_col[myId].size()) {
                    //                    int d_row = res_packages_same_col[myId][ind].d_row;
                    int d_row = ROW(res_packages_same_col[myId][ind].dst_id);
                    int s = 1;
                    counter[d_row].push_back(s);
                }
            }
            nFirstPr = counter.size();

            putc_schedules_same_col[myId].push_back(nFirstPr);
        }


        ///////generate getc_schedules
        for (int myId = 0; myId < nThread; myId++) {

            //count how many getrputc will do in my thread 
            int myRow = ROW(myId);
            int myCol = COL(myId);
            getc_schedules_same_col[myId].push_back(0);

            //loop the threads in my column, find which will putc to myId
            for (int row = 0; row < 8; row++) {
                int mySameColId = myCol + row * 8;
                if (mySameColId != myId) {//loop the threads in my column 
                    int nFirstPr = putc_schedules_same_col[mySameColId][cycle];
                    for (int n = 0; n < nFirstPr; n++) {
                        int ind = nextCycleStartPos[mySameColId] + n;
                        //                        if (res_packages_same_col[mySameColId][ind].d_row == myRow)
                        if (ROW(res_packages_same_col[mySameColId][ind].dst_id) == myRow) {
                            getc_schedules_same_col[myId][cycle]++;
                            //                            getrputc_list[myId].push_back(res_packages[mySameColId][ind].d_row); //record which row myId will putc
                        }
                    }
                }
            }
        }
        //
        //
        for (int myId = 0; myId < nThread; myId++) {
            int nFirstPr = putc_schedules_same_col[myId][cycle];
            nextCycleStartPos[myId] += nFirstPr;
        }

        cycle++;
    }

    //    cout << "\nGenerate column schedules done!" << endl;
    //FOR TEST///////////
    for (int myId = 0; myId < nThread; myId++) {

        if (putc_schedules_same_col[myId].size() != getc_schedules_same_col[myId].size()) {
            cout << "error" << endl;
            abort();
        }
    }
}

void RlmpiInitializer::generate_data_same_col() {
    srand(time(NULL));
    for (int myId = 0; myId < nThread; myId++) {

        int myCol = COL(myId);
        int myRow = ROW(myId);

        for (int myRow = 0; myRow < 8; myRow++) {
            int dstId = myCol + myRow * 8; //the id in my column

            if (dstId != myId) {
                int nSend = rand() % maxNPack;
                //                int nSend = 1;
                for (int i = 0; i < nSend; i++) {
                    Pack pack;
                    pack.src_id = myId;
                    pack.dst_id = dstId;
                    //                    pack.d_col = COL(pack.dst_id);
                    //                    pack.d_row = ROW(pack.dst_id);
                    //                    pack.s_col = COL(pack.src_id);
                    //                    pack.s_row = ROW(pack.src_id);
                    pack.data[0] = double(myId) + double(i) / 1e2;
                    pack.data[1] = double(myId) + double(i) / 1e3;
                    pack.data[2] = double(myId) + double(i) / 1e4;
                    packages[myId][dstId].push_back(pack);
                }
            }
        }
    }
    //    cout << "Packages generated" << endl;
}

void RlmpiInitializer::reorder_packages_same_col() {

    for (int myId = 0; myId < nThread; myId++) {

        vector<Pack> all_packs; // All packages sent of myId

        int sThread[nThread]; // counter of packages of dstId
        for (int i = 0; i < nThread; i++) {
            sThread[i] = 0;
        }
        for (int n = 0; n < table[myId].nPutcSameCol; n++) {
            for (map<int, vector<Pack> >::iterator it = same_col_packages[myId].begin();
                    it != same_col_packages[myId].end(); it++) {

                int dstId = it->first;
                int nSend = it->second.size(); // n packages sent to dstId
                const vector<Pack> &packs = it->second; //packages to dstId

                if (sThread[dstId] < nSend) {
                    all_packs.push_back(packs[sThread[dstId]++]); //see the excel in folder
                }
            }
        }
        res_packages_same_col[myId] = all_packs;
    }

    //    cout << "same column package reorder done!" << endl;
}

void RlmpiInitializer::generate_table_same_col() {
    for (int myId = 0; myId < nThread; myId++) {
        table[myId].nGetcSameCol = 0;
        table[myId].nPutcSameCol = 0;
    }

    int dst_list[64][64];
    int sendN_list[64][64];
    for (int myId = 0; myId < 64; myId++) {
        for (int dstId = 0; dstId < 64; dstId++) {
            sendN_list[myId][dstId] = 0;
            dst_list[myId][dstId] = -1;
        }
    }

    for (int myId = 0; myId < nThread; myId++) {

        for (map<int, vector<Pack> >::iterator it = same_col_packages[myId].begin();
                it != same_col_packages[myId].end(); it++) {

            const vector<Pack> &packs = it->second;
            int dstId = it->first;
            int nSend = it->second.size();

            dst_list[myId][dstId] = dstId;
            sendN_list[myId][dstId] = nSend;
        }
    }

    generate_register_transform_table_same_col(dst_list, sendN_list, table);

    //    cout << "Generate same column table done" << endl;
}

void RlmpiInitializer::generate_dst_sequence() {

    for (int myId = 0; myId < nThread; myId++) {
        dst_sequence[myId][0] = myId;
        int s = 1;
        for (int i = 1; i < nThread; i++) {
            int upperId = myId + i;
            int downId = myId - i;

            if (upperId <= 63 && downId >= 0) {
                dst_sequence[myId][s] = upperId;
                s++;
                dst_sequence[myId][s] = downId;
                s++;
            }
            if (upperId <= 63 && downId < 0) {
                dst_sequence[myId][s] = upperId;
                s++;
            }
            if (upperId > 63 && downId >= 0) {

                dst_sequence[myId][s] = downId;
                s++;
            }
            //            cout << i << endl;

        }
    }
    ///FOR TEST
    //    for (int myId = 0; myId < nThread; myId++) {
    //        for (int i = 0; i < nThread; i++) {
    //            cout << "myId " << myId << " " << dst_sequence[myId][i] << endl;
    //        }
    //    }
}

void RlmpiInitializer::generate_data_skew(const vector<vector <int8LDM> > &src2dst, const vector<vector <int8LDM> > &nData) {
    for (int myId = 0; myId < nThread; myId++) {
        packages[myId].clear();
        for (int s = 0; s < src2dst[myId].size(); s++) {
            Pack pack;
            pack.src_id = myId;
            int dstId = src2dst[myId][s];
            pack.dst_id = dstId;

            pack.cva = nData[myId][s];

            pack.data[0] = 0;
            pack.data[1] = 0;
            pack.data[2] = 0;

            packages[myId][dstId].push_back(pack);
        }
    }
    //    cout << "Not on same row column packages generated" << endl;
}

void RlmpiInitializer::assemble_packages() {
    all_res_packages.resize(nThread);

    for (int myId = 0; myId < nThread; myId++) {
        //        DISP2(res_packages[myId].size(), res_packages_same_row[myId].size());
        if (res_packages[myId].size() > 0)
            all_res_packages[myId].insert(all_res_packages[myId].end(), res_packages[myId].begin(), res_packages[myId].end());

        if (res_packages_same_row[myId].size() > 0)
            all_res_packages[myId].insert(all_res_packages[myId].end(), res_packages_same_row[myId].begin(), res_packages_same_row[myId].end());

        if (res_packages_same_col[myId].size() > 0) {
            all_res_packages[myId].insert(all_res_packages[myId].end(), res_packages_same_col[myId].begin(), res_packages_same_col[myId].end());
        }
    }
    generate_recv_position();
}

void RlmpiInitializer::generate_recv_position() {

    for (int myId = 0; myId < nThread; myId++) {

        int s = 0;
        map<int, vector<int> > qq;
        for (int srcId = 0; srcId < nThread; srcId++) {
            for (int n = 0; n < all_res_packages[srcId].size(); n++) {
                int dst_id = all_res_packages[srcId][n].dst_id;
                int src_id = all_res_packages[srcId][n].src_id;
                if (myId == dst_id && src_id != myId) {
                    qq[src_id].push_back(s);
                    s++;
                }
            }
        }

        //setup the position
        for (map<int, vector<int> >::iterator it = qq.begin(); it != qq.end(); it++) {
            int src_id = it->first;
            int s = 0;
            for (int m = 0; m < all_res_packages[src_id].size(); m++) {
                if (all_res_packages[src_id][m].dst_id == myId) {
                    all_res_packages[src_id][m].res_pos = it->second[s];
                    s++;
                    if (myId == 63) {
                        //                        DISP2(s,all_res_packages[0].size());
                    }
                }
            }
        }

    }

}

void RlmpiInitializer::init(const vector<vector<int32> > &sendDstLists) {

    for (int myId = 0; myId < nThread; myId++) {
        vector<int> res = get_destination_pool(myId);
        destination_pool.push_back(res);
    }
    packages.resize(nThread);
    same_col_packages.resize(nThread);
    same_row_packages.resize(nThread);

    res_packages.resize(nThread);
    res_packages_same_row.resize(nThread);
    res_packages_same_col.resize(nThread);

    putc_schedules_same_col.resize(nThread);
    getc_schedules_same_col.resize(nThread);

    putr_schedules_same_row.resize(nThread);
    getr_schedules_same_row.resize(nThread);

    generate_dst_sequence();

    same_row_dst.resize(nThread);
    same_col_dst.resize(nThread);
    not_col_row_dst.resize(nThread);

    for (int myId = 0; myId < nThread; myId++) {
        for (int n = 0; n < sendDstLists[myId].size(); n++) {
            int dst = sendDstLists[myId][n];
            int d_col = COL(dst);
            int d_row = ROW(dst);
            int my_col = COL(myId);
            int my_row = ROW(myId);
            if (dst == myId) {
                DISP(" Can't Send to oneself");
                DISP2(myId, n);
                abort();
            }
            if (d_col != my_col && d_row != my_row) {
                not_col_row_dst[myId].push_back(dst);
            }
            if (d_col == my_col && d_row != my_row) {
                same_col_dst[myId].push_back(dst);
            }
            if (d_row == my_row && d_col != my_col) {
                same_row_dst[myId].push_back(dst);
            }
        }
    }

    not_col_row_Ndata.resize(nThread);
    same_col_Ndata.resize(nThread);
    same_row_Ndata.resize(nThread);

    for (int myId = 0; myId < nThread; myId++) {

        for (int n = 0; n < same_col_dst[myId].size(); n++) {
            same_col_Ndata[myId].push_back(6);
        }

        for (int n = 0; n < same_row_dst[myId].size(); n++) {
            same_row_Ndata[myId].push_back(6);
        }

        for (int n = 0; n < not_col_row_dst[myId].size(); n++) {
            not_col_row_Ndata[myId].push_back(6);
        }
    }
    generate_data_skew(not_col_row_dst, not_col_row_Ndata);
    generate_table_skew();
    reorder_packages_skew();
    generate_schedule_skew();

    generate_data_same_row(same_row_dst, same_row_Ndata);
    generate_table_same_row();
    reorder_packages_same_row();
    generate_schedule_same_row();

    generate_data_same_col(same_col_dst, same_col_Ndata);
    generate_table_same_col();
    reorder_packages_same_col();
    generate_schedule_same_col();
    assemble_packages();
}

void RlmpiInitializer::init(const vector<vector<int32> > &sendDstLists, const vector<vector<int32> > &nDataLists) {
    for (int myId = 0; myId < nThread; myId++) {
        vector<int> res = get_destination_pool(myId);
        destination_pool.push_back(res);
    }
    packages.resize(nThread);
    same_col_packages.resize(nThread);
    same_row_packages.resize(nThread);

    res_packages.resize(nThread);
    res_packages_same_row.resize(nThread);
    res_packages_same_col.resize(nThread);

    putc_schedules_same_col.resize(nThread);
    getc_schedules_same_col.resize(nThread);

    putr_schedules_same_row.resize(nThread);
    getr_schedules_same_row.resize(nThread);

    generate_dst_sequence();

    same_row_dst.resize(nThread);
    same_col_dst.resize(nThread);
    not_col_row_dst.resize(nThread);

    not_col_row_Ndata.resize(nThread);
    same_col_Ndata.resize(nThread);
    same_row_Ndata.resize(nThread);


    for (int myId = 0; myId < nThread; myId++) {
        for (int n = 0; n < sendDstLists[myId].size(); n++) {
            int dst = sendDstLists[myId][n];
            int d_col = COL(dst);
            int d_row = ROW(dst);
            int my_col = COL(myId);
            int my_row = ROW(myId);
            if (dst == myId) {
                DISP(" Can't Send to oneself");
                DISP2(myId, n);
                abort();
            }
            if (d_col != my_col && d_row != my_row) {
                not_col_row_dst[myId].push_back(dst);
                not_col_row_Ndata[myId].push_back(nDataLists[myId][n]);

            }
            if (d_col == my_col && d_row != my_row) {
                same_col_dst[myId].push_back(dst);
                same_col_Ndata[myId].push_back(nDataLists[myId][n]);

            }
            if (d_row == my_row && d_col != my_col) {
                same_row_dst[myId].push_back(dst);
                same_row_Ndata[myId].push_back(nDataLists[myId][n]);
            }
        }
    }



    generate_data_skew(not_col_row_dst, not_col_row_Ndata);
    generate_table_skew();
    reorder_packages_skew();
    generate_schedule_skew();

    generate_data_same_row(same_row_dst, same_row_Ndata);
    generate_table_same_row();
    reorder_packages_same_row();
    generate_schedule_same_row();

    generate_data_same_col(same_col_dst, same_col_Ndata);
    generate_table_same_col();
    reorder_packages_same_col();
    generate_schedule_same_col();
    assemble_packages();
}

vector<int> RlmpiInitializer::GetLDM_MemoryUsage() {
    vector<int> mem_len(nThread, 0);
    for (int myId = 0; myId < nThread; myId++) {
        mem_len[myId] += all_res_packages[myId].size() * sizeof (Pack)*2;
        mem_len[myId] += sizeof (Table);
        mem_len[myId] += sizeof (int8LDM) * getc_schedules_same_col[myId].size()*2;
        mem_len[myId] += sizeof (int8LDM) * getr_schedules_same_row[myId].size()*2;
        mem_len[myId] += sizeof (int8LDM) * getc_schedules[myId].size()*3;
    }
    return mem_len;
}

void RlmpiInitializer::copyRlmpiInfo(RlmpiInfo *schedule_data) {
    schedule_data->nCycle = putr_schedules[0].size();
    schedule_data->nCycleSameRow = putr_schedules_same_row[0].size();
    schedule_data->nCycleSameCol = putc_schedules_same_col[0].size();


    if (schedule_data->nCycle > MaxNCycle) {
        DISP2("nCycle ", schedule_data->nCycle);
        abort();
    }
    if (schedule_data->nCycleSameRow > MaxNCycle) {
        DISP2("nCycleSameRow ", schedule_data->nCycleSameRow);
        abort();
    }
    if (schedule_data->nCycleSameCol > MaxNCycle) {
        DISP2("nCycleSameCol ", schedule_data->nCycleSameCol);
        abort();
    }

    for (int myId = 0; myId < nThread; myId++) {

        schedule_data->table[myId] = table[myId];
        int total_send = table[myId].nPutrSkew + table[myId].nPutcSameCol + table[myId].nPutrSameRow;
        int total_recv = table[myId].nGetcSkew + table[myId].nGetcSameCol + table[myId].nGetrSameRow;


        if (total_send > MaxNPackages) {
            DISP(total_send);
            abort();
        }




        schedule_data->package[myId] = new Pack[total_send];
        schedule_data->dstId_list[myId] = new int8LDM[total_send];
        schedule_data->srcId_list[myId] = new int8LDM[total_send];
        schedule_data->resPos_list[myId] = new int16LDM[total_send];
        schedule_data->cva_list[myId] = new int16LDM[total_send];
        schedule_data->cvb_list[myId] = new int16LDM[total_send];
        //init header by copy packages
        for (int j = 0; j < total_send; j++) {
            schedule_data->package[myId][j] = all_res_packages[myId][j];

            schedule_data->dstId_list[myId][j] = all_res_packages[myId][j].dst_id;
            schedule_data->srcId_list[myId][j] = all_res_packages[myId][j].src_id;
            schedule_data->resPos_list[myId][j] = all_res_packages[myId][j].res_pos;
            schedule_data->cva_list[myId][j]=all_res_packages[myId][j].cva;
            schedule_data->cvb_list[myId][j]=all_res_packages[myId][j].cvb;
        }

        ///////load schedules
        schedule_data->getc_schedules[myId] = new int8LDM[getc_schedules[myId].size()];
        schedule_data->getrputc_schedules[myId] = new int8LDM[getrputc_schedules[myId].size()];
        schedule_data->putr_schedules[myId] = new int8LDM[putr_schedules[myId].size()];

        for (int j = 0; j < getc_schedules[myId].size(); j++) {
            schedule_data->getc_schedules[myId][j] = getc_schedules[myId][j];
            schedule_data->getrputc_schedules[myId][j] = getrputc_schedules[myId][j];
            schedule_data->putr_schedules[myId][j] = putr_schedules[myId][j];
        }

        schedule_data->getr_schedules_same_row[myId] = new int8LDM[getr_schedules_same_row[myId].size()];
        schedule_data->putr_schedules_same_row[myId] = new int8LDM[putr_schedules_same_row[myId].size()];
        for (int j = 0; j < putr_schedules_same_row[myId].size(); j++) {
            schedule_data->getr_schedules_same_row[myId][j] = getr_schedules_same_row[myId][j];
            schedule_data->putr_schedules_same_row[myId][j] = putr_schedules_same_row[myId][j];
        }

        schedule_data->getc_schedules_same_col[myId] = new int8LDM[getc_schedules_same_col[myId].size()];
        schedule_data->putc_schedules_same_col[myId] = new int8LDM[putc_schedules_same_col[myId].size()];
        for (int j = 0; j < putc_schedules_same_col[myId].size(); j++) {
            schedule_data->getc_schedules_same_col[myId][j] = getc_schedules_same_col[myId][j];
            schedule_data->putc_schedules_same_col[myId][j] = putc_schedules_same_col[myId][j];
        }
    }



}
#endif
