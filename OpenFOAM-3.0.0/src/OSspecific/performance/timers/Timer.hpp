#ifndef TIMER_HPP
#define TIMER_HPP

#include <map>
#include <string>

class swTimer
{
public:

    static std::map<std::string, double> timeSum;
    static std::map<std::string, double> timeStart;
    static std::map<std::string, double> timeEnd;
    static std::map<std::string, int>    count;

    void static startTimer(std::string in);
    static void endTimer(std::string in);

    ~swTimer();

    static void printTimer();
    static void printTimer(std::string in);

    static void maxTimeSum();
};

#if(SWTIMER)
    #define TIMER_START(name) \
    swTimer::startTimer(name);
#else
    #define TIMER_START(name)
#endif

#if(SWTIMER)
    #define TIMER_END(name) \
    swTimer::endTimer(name);
#else
    #define TIMER_END(name)
#endif

#if(SWTIMER)
    #define TIMER_PRINT() \
    swTimer::printTimer();
#else
    #define TIMER_PRINT()
#endif


#endif
