#ifndef SUNWAYMACROS_H
#define SUNWAYMACROS_H

#include "mpi.h"

typedef long swInt64;
typedef int swInt32;
typedef int swInt;
typedef double swFloat64;
typedef float swFloat32;
typedef double swFloat;

#define MAX_CARED_RERR (1e-9)
#define MIN_FP (1e-15)
// this macro takes two arrays, in which the element
// can be the operand of "mag" that returns a scalar.
#define HOST2SLAVE_MPI_CHECK(var1_, var2_) \
{ \
    double max_rerr=0; \
    forAll(var1_, celli) \
    { \
        double mag1 = mag(var1_[celli]); \
        double mag2 = mag(var2_[celli]); \
        if(mag1 != 0) \
        { \
            double rerr=fabs((mag2-mag1)/mag1); \
            max_rerr = (max_rerr > rerr) ? max_rerr : rerr; \
        } \
        else if(mag2 > MIN_FP) \
        { \
            max_rerr = 1; \
            break; \
        } \
    } \
    double max_rerr_buf = max_rerr; \
    MPI_Reduce(&max_rerr_buf, \
               &max_rerr, \
               1, \
               MPI_DOUBLE, \
               MPI_SUM, \
               0, \
               MPI_COMM_WORLD \
              ); \
    int rank; \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
    if(rank == 0) \
    if(max_rerr > MAX_CARED_RERR) \
        std::cout<<"***Error at " \
                 <<__FILE__<<__LINE__<<":"<<std::endl \
                 <<"max relative error:" \
                 <<max_rerr<<", exceeds "<<MAX_CARED_RERR \
                 <<std::endl; \
}

// this macro takes two arrays, in which the element
// can be the operand of "mag" that returns a scalar.
#define HOST2SLAVE_CHECK(var1_, var2_) \
{ \
    double max_rerr=0; \
    forAll(var1_, celli) \
    { \
        double mag1 = mag(var1_[celli]); \
        double mag2 = mag(var2_[celli]); \
        if(mag1 != 0) \
        { \
            double rerr=fabs((mag2-mag1)/mag1); \
            max_rerr = (max_rerr > rerr) ? max_rerr : rerr; \
        } \
        else if(mag2 > MIN_FP) \
        { \
            max_rerr = 1; \
            break; \
        } \
    } \
    if(max_rerr > MAX_CARED_RERR) \
        std::cout<<"***Error at " \
                 <<__FILE__<<__LINE__<<":"<<std::endl \
                 <<"max relative error:" \
                 <<max_rerr<<", exceeds "<<MAX_CARED_RERR \
                 <<std::endl; \
}

// this macro takes two arrays, in which the element
// can be the operand of "mag" that returns a scalar.
#define HOST2SLAVE_CHECK_DETAIL(var1_, var2_) \
{ \
    double max_rerr=0; \
    forAll(var1_, celli) \
    { \
        double mag1 = mag(var1_[celli]); \
        double mag2 = mag(var2_[celli]); \
        if(mag1 != 0) \
        { \
            double rerr=fabs((mag2-mag1)/mag1); \
            max_rerr = (max_rerr > rerr) ? max_rerr : rerr; \
            if(rerr > MAX_CARED_RERR) \
                std::cout<<"mag1="<<mag1 \
                         <<", mag2="<<mag2<<std::endl; \
        } \
        else if(mag2 > MIN_FP) \
        { \
            std::cout<<"mag1="<<mag1 \
                     <<", mag2="<<mag2<<std::endl; \
            max_rerr = 1; \
            break; \
        } \
    } \
    if(max_rerr > MAX_CARED_RERR) \
        std::cout<<"***Error at " \
                 <<__FILE__<<__LINE__<<":"<<std::endl \
                 <<"max relative error:" \
                 <<max_rerr<<", exceeds "<<MAX_CARED_RERR \
                 <<std::endl; \
}

// this macro takes two scalars
#define HOST2SLAVE_MPI_COMPARE(mag1, mag2) \
{ \
    double rerr = 0; \
    if(mag1 != 0) \
    { \
        rerr=fabs((mag2-mag1)/mag1); \
    } \
    else if(fabs(mag2) > MIN_FP) \
    { \
        rerr = 1.0; \
    } \
    double rerr_buf = rerr; \
    MPI_Reduce(&rerr_buf, \
               &rerr, \
               1, \
               MPI_DOUBLE, \
               MPI_SUM, \
               0, \
               MPI_COMM_WORLD \
              ); \
    int rank; \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
    if(rank == 0) \
    if(rerr > MAX_CARED_RERR) \
        std::cout<<"***Error at " \
                 <<__FILE__<<__LINE__<<":"<<std::endl \
                 <<"relative error:" \
                 <<rerr<<", exceeds "<<MAX_CARED_RERR \
                 <<std::endl; \
}

// this macro takes two scalars
#define HOST2SLAVE_COMPARE(mag1, mag2) \
{ \
    double rerr = 0; \
    if(mag1 != 0) \
    { \
        rerr=fabs((mag2-mag1)/mag1); \
        if(rerr > MAX_CARED_RERR) \
            std::cout<<"mag1="<<mag1 \
                     <<", mag2="<<mag2<<std::endl; \
    } \
    else if(fabs(mag2) > MIN_FP) \
    { \
        rerr = 1.0; \
        std::cout<<"mag1="<<mag1 \
                 <<", mag2="<<mag2<<std::endl; \
    } \
    if(rerr > MAX_CARED_RERR) \
        std::cout<<"***Error at " \
                 <<__FILE__<<__LINE__<<":"<<std::endl \
                 <<"relative error:" \
                 <<rerr<<", exceeds "<<MAX_CARED_RERR \
                 <<std::endl; \
}

#endif
