#ifndef _COMMON_H
#define _COMMON_H

#include "sort.h"

#ifdef DEBUG
#define _DEBUG_DEBUG 1
#else
#define _DEBUG_DEBUG 0
#endif

#define Dprintf(fmt, ...) \
        do { if (_DEBUG_DEBUG) \
                fprintf(stderr, "%s:%d: %s(): " fmt, \
                    __FILE__, \
                    __LINE__, \
                    __func__, \
                    ##__VA_ARGS__);\
        } while (0)

#define NDprintf(...) \
        do { if (!_DEBUG_DEBUG) \
                fprintf(stderr, ##__VA_ARGS__);\
        } while (0)

#define taskfor(j, n, ...) for(int _i = 0; _i < n; ) \
{ \
  int lim = _i+(n)/num_procs; \
  if(lim < _i+1) \
    lim = _i+1; \
  if(lim > (n)) \
    lim = (n); \
  j = _i; \
  PRAGMA(omp task if((n) > for_ser_n) firstprivate(j, lim, ##__VA_ARGS__) untied) \
  { \
    for(; j < lim; ++j)

#define endtaskfor } _i=lim;} \
  PRAGMA(omp taskwait)

#define PRAGMA(x) _Pragma(#x)

long long randull(unsigned int *seed);
void psum(int *data, int n, int *data2);

#endif // _COMMON_H
