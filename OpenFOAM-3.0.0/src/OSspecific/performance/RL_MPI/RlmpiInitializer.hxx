/* 
 * File:   RegisterLevelMsgPass.hxx
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
#include "RlmpiSharedType.h"
#ifdef SUNWAY
extern "C" {
#include <athread.h>
}
#endif
using namespace std;
#ifndef DISP
#define DISP(x) std::cout << __FILE__<<" "<<__LINE__<< ",    " << #x ": " << (x) << std::endl
#endif


#ifndef DISP2
#define DISP2(x,y) std::cout << __FILE__<<" "<<__LINE__<< ",    " << #x ": " << x <<" "<<y<< std::endl
#endif
static inline double timestamp() {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return ts.tv_sec + ts.tv_nsec * 1e-9;
}

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
        void generate_data();   // first step is generate data
        void generate_data_un_col_row(const vector<vector < int8LDM> > &not_col_row_dst); // first step is generate data
        void generate_data_un_col_row(const vector<vector <int8LDM> > &not_col_row_dst, const vector<vector <int8LDM> > &not_col_row_Ndata); // first step is generate data
        void generate_table();  // generate table for every thread
        void reorder_packages(); // reorder the sequential of packages
        void generate_schedule(); 
        void reorder_packages2();

        void generate_data_same_row(); // first step is generate data
        void generate_data_same_row(const vector<vector <int8LDM> > &not_col_row_dst, const vector<vector <int8LDM> > &nData);
        void reorder_packages_same_row();
        void generate_table_same_row();
        void generate_schedule_same_row();

        ////////
        void generate_data_same_col(); // first step is generate data
        void generate_data_same_col(const vector<vector <int8LDM> > &not_col_row_dst, const vector<vector <int8LDM> > &nData);
        void reorder_packages_same_col();
        void generate_table_same_col();
        void generate_schedule_same_col();
        ///////
        void assemble_packages();


        ///////////////////////
        friend void generate_register_transform_table(int(*dst_list)[64],
                int(*sendN)[64], Table *table);
        friend void generate_register_transform_table_same_row(const int(*dst_list)[64],
                const int(*sendN)[64], Table *table);
        friend void generate_register_transform_table_same_col(int(*dst_list)[64],
                int(*sendN)[64], Table *table);
        void transpose_matrix(vector<Pack> &pack);
        void write_packages();
        void write_schedule();

        void generate_recv_position();
        void init(const vector<vector<ThreadID> > &sendDstLists);
        void copyinfo(Schedule *reg_data);
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




