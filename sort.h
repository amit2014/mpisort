#ifndef _SORT_H
#define _SORT_H

#include "payloadsize.h"  // Defines LOADSIZE
typedef enum {
   QUICK,
   MERGE,
   RADIX,
   BEST
}  SortType;

typedef struct {
   long long *key;
   char payload[LOADSIZE];
} dataType;

void pSort(dataType *data, int ndata, SortType sorter=BEST);

#endif // _SORT_H