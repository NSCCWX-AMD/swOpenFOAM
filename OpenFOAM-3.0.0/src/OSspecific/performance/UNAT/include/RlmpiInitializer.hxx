/* 
 * File:   RlmpiInitializer.hxx
 * Author: Geng Chen (gengchn@gmail.com)
 *
 * Created on November 30, 2017, 3:49 PM
 */

#ifndef FILE_RLMPI_INITIALLIZER_H
#define FILE_RLMPI_INITIALLIZER_H
#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)

#include <vector>
#include <map>
#include <time.h>
// #include <iostream>
#include "RlmpiShared.h"
#ifdef SUNWAY
extern "C" {
#include <athread.h>
}
#endif
using namespace std;
#ifndef DISP
#define DISP(x)
//#define DISP(x) std::cout << __FILE__<<" "<<__LINE__<< ",    " << #x ": " << (x) << std::endl
#endif


#ifndef DISP2
#define DISP2(x,y)
//#define DISP2(x,y) std::cout << __FILE__<<" "<<__LINE__<< ",    " << #x ": " << x <<" "<<y<< std::endl
#endif

static inline double timestamp() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// template < class T >
// std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
//     //    os << "[";
//     for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end() - 1; ++it) {
//         os << "  " << *it << ", ";
//     }
//     os << *v.rbegin();
//     os << ";";
//     return os;
// }

class RlmpiInitializer {
public:
    Table table[64]; //length is 64
    static const int nThread = 64;
    static const int maxNPack = 10;
    static const int maxNdst = 49;

    static const int nSameRow = 7;
    static const int nSameCol = 7;

    vector<vector<int8LDM> > putr_schedules;
    vector<vector<int8LDM> > getrputc_schedules;
    vector<vector<int8LDM> > getc_schedules;
    vector<vector<Pack> > res_packages;
    vector<vector<Pack> > res_packages_same_col;
    vector<vector<Pack> > res_packages_same_row;
    vector<vector<Pack> > all_res_packages;
    ///////////
    vector<vector<int8LDM> > putr_schedules_same_row;
    vector<vector<int8LDM> > getr_schedules_same_row;

    ///////////
    vector<vector<int8LDM> > putc_schedules_same_col;
    vector<vector<int8LDM> > getc_schedules_same_col;

public:
    RlmpiInitializer(); //constructor
    void init(const vector<vector<int32> > &sendDstLists);
    void init(const vector<vector<int32> > &sendDstLists,const vector<vector<int32> > &nDataLists);
    void set_cva(const vector<vector<int32> > &cvaLists);
    void set_cvb(const vector<vector<int32> > &cvbLists);
    void get_package_len(int myId);
    void copyRlmpiInfo(RlmpiInfo *rlmpi_info);
    vector<int> GetLDM_MemoryUsage();
    void write_packages();
    void write_schedule();

protected:
    void generate_data(); // first step is generate data
    void generate_data_skew(const vector<vector <int8LDM> > &not_col_row_dst); // first step is generate data
    void generate_data_skew(const vector<vector <int8LDM> > &not_col_row_dst,
            const vector<vector <int8LDM> > &not_col_row_Ndata); // first step is generate data
    void generate_table_skew(); // generate table for every thread
    void reorder_packages_skew(); // reorder the sequential of packages
    void generate_schedule_skew();
    void reorder_packages_skew2();

    void generate_data_same_row(); // first step is generate data
    void generate_data_same_row(const vector<vector <int8LDM> > &not_col_row_dst,
            const vector<vector <int8LDM> > &nData);
    void reorder_packages_same_row();
    void generate_table_same_row();
    void generate_schedule_same_row();

    ////////
    void generate_data_same_col(); // first step is generate data
    void generate_data_same_col(const vector<vector <int8LDM> > &not_col_row_dst,
            const vector<vector <int8LDM> > &nData);
    void reorder_packages_same_col();
    void generate_table_same_col();
    void generate_schedule_same_col();
    ///////
    void assemble_packages();


    ///////////////////////
    friend void generate_register_transform_table_skew(int(*dst_list)[64],
            int(*sendN)[64], Table *table);
    friend void generate_register_transform_table_same_row(const int(*dst_list)[64],
            const int(*sendN)[64], Table *table);
    friend void generate_register_transform_table_same_col(int(*dst_list)[64],
            int(*sendN)[64], Table *table);

    void generate_recv_position();
    ///////////////////////

protected:
    vector<int> get_destination_pool(int myId);
    void generate_dst_sequence();
    vector<vector<int> > destination_pool;
    vector<map<int, vector<Pack> > > packages;
    vector<map<int, vector<Pack> > > non_same_col_row_packages;
    vector<map<int, vector<Pack> > > same_col_packages;
    vector<map<int, vector<Pack> > > same_row_packages;
    vector<map<int, vector<Pack> > > all_packages;

    int dst_sequence[64][64];
    static const int nRecReg = 7;
    static const int nSendReg = 7;


    //define how many useful data in each package
    vector<vector<int8LDM> > same_row_Ndata;
    vector<vector<int8LDM> > same_col_Ndata;
    vector<vector<int8LDM> > not_col_row_Ndata;

    vector<vector<int8LDM> > same_row_dst;
    vector<vector<int8LDM> > same_col_dst;
    vector<vector<int8LDM> > not_col_row_dst;
};

#endif /* REGISTERLEVELMSGPASS_H */




