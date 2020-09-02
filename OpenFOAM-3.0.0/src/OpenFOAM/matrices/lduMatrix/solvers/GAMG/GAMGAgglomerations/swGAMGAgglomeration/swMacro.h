#ifndef SWMACRO_H
#define SWMACRO_H

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCKNUM64K 64
#define EPS 1e-6
#define DEBUG

#ifdef DEBUG
#define LOG(format,...) printf("File: "__FILE__",Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)
#else
#define LOG(format,...)
#endif

#define SLAVE_FUNC(funcname) slave_##funcname
#define ALIGNED(addr) ((( ( (unsigned long)(addr) - 1)>>5)+1)<<5)
#define ArraySize 57344

typedef int    swInt;
typedef int    swInt32;
typedef long   swInt64;

typedef double swFloat;
typedef float  swFloat32;
typedef double swFloat64;

#ifdef __cplusplus
}
#endif

#endif
