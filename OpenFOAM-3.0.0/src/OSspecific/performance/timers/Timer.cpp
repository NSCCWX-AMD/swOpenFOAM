#include <iostream>
#include <sys/time.h>
#include "mpi.h"
#include "Timer.hpp"
#include "Timer_c.h"

using std::map;
using std::string;
using std::cout;
using std::endl;

map<string, double> swTimer::timeSum;
map<string, double> swTimer::timeStart;
map<string, double> swTimer::timeEnd;
map<string, int>    swTimer::count;

inline double getSystemTime()
{
    struct timeval timer;
    gettimeofday(&timer, 0);
    return ((double)(timer.tv_sec) + (double)(timer.tv_usec)*1.0e-6);
}

void swTimer::startTimer(string in)
{
    timeStart[in] = getSystemTime();
}


void swTimer::endTimer(string in)
{
    timeEnd[in] = getSystemTime();
    map<string, double>::iterator it = timeSum.find(in);
    if(it==timeSum.end())
    {
        timeSum[in] = 0;
        count[in] = 0;
    }
    timeSum[in] += timeEnd[in] - timeStart[in];
    count[in]++;
}


void swTimer::printTimer(string in)
{
    map<string, double>::iterator it = timeSum.find(in);

    if(it==timeEnd.end())
    {
        cout << "there is no swTimer. \"" << in << endl;
    }
    else
    {
        cout << it->first <<" Time: "
             << it->second <<"s, count: "
             << count[it->first] << ",  average: "
             << it->second/count[it->first] << "s" << endl;
    }
}

void swTimer::printTimer()
{
    int initialized;
    MPI_Initialized(&initialized);

    int rootProc = 0;

    if(initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rootProc);
    }

    if(!rootProc)
    {
        map<string, double>::iterator it;

        for(it=timeSum.begin(); it!=timeSum.end(); ++it)
        {
            cout << it->first << " Time: "
                 << it->second << "s, count: "
                 << count[it->first] << ",  average: "
                 << it->second/count[it->first] << "s" << endl;
        }
    }
}


void swTimer::maxTimeSum()
{
    map<string, double>::iterator it;

    for(it=timeSum.begin(); it!=timeSum.end(); ++it)
    {
        double time1 = it->second, time2;
        MPI_Allreduce(&time1, &time2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        it->second = time2;
    }
}


swTimer::~swTimer()
{
    map<string, double>::iterator it;
    std::cout << "swTimer ~~~" << endl;

    for(it=timeSum.begin(); it!=timeSum.end(); ++it)
    {
        cout << it->first << " Time: "
             << it->second << "s, count: "
             << count[it->first] << ",  average: "
             << it->second/count[it->first] << "s" << endl;
    }
}


//- interfaces for c
void swTimerStart(const char* name)
{
    swTimer::startTimer(name);
}


void swTimerEnd(const char* name)
{
    swTimer::endTimer(name);
}


void swTimerPrintAll()
{
    swTimer::printTimer();
}


void swTimerPrint(const char* name)
{
    swTimer::printTimer(name);
}


void swTimerMaxSum()
{
    swTimer::maxTimeSum();
}



//- interfaces for Fortran
void swtimerstart_(const char* name)
{
    swTimer::startTimer(name);
}


void swtimerend_(const char* name)
{
    swTimer::endTimer(name);
}


void swtimerprintall_()
{
    swTimer::printTimer();
}


void swtimerprint_(const char* name)
{
    swTimer::printTimer(name);
}


void swtimermaxsum_()
{
    swTimer::maxTimeSum();
}
