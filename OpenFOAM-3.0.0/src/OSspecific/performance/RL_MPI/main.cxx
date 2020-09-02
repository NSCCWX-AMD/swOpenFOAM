#include "RlmpiInitializer.hxx"
#include "iostream"
#include <vector>
using namespace std;
#ifdef SUNWAY
extern "C" {
#include <athread.h>
#include <algorithm>
    void slave_initTable(Schedule*data);
    void slave_test(Schedule*data);
    void slave_test_athread_get(double* in);
}
#endif

int main() {
#ifdef SUNWAY
    athread_init();
#endif
	Schedule *schedule_data = new Schedule[40];
	for(int idx=0;idx<40;idx++){
		RlmpiInitializer reg;
	    vector<vector<ThreadID> > sendList(64);
		for (int i = 0; i < 64; i++) {
			int ntimes = 2;
	        //        sendList[i].resize(63 * ntimes);
		    int s = 0;
			for (int b = 0; b < ntimes; b++) {
				for (int j = 0; j < 64; j++) {
					if (j != i) {
						sendList[i].push_back(j);
	                }
		        }
			}
	    }
	    reg.init(sendList);
		reg.copyinfo(&schedule_data[idx]);
	    schedule_data[idx].destroy = 0;
	}

#ifdef SUNWAY
	for(int idx=0;idx<40;idx++){
		printf("%d\n",idx);
	    schedule_data[idx].destroy = 0;
		__real_athread_spawn((void *) slave_initTable, &schedule_data[idx]);
	    athread_join();
#endif

#ifdef SUNWAY
    __real_athread_spawn((void *) slave_test, &schedule_data[idx]);
    athread_join();
#endif
	    schedule_data[idx].destroy = 1;
		__real_athread_spawn((void *) slave_initTable, &schedule_data[idx]);
	    athread_join();
	}
 



#ifdef SUNWAY
    athread_halt();
#endif
}
