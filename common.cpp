#include <stdlib.h>
#include <omp.h>
#include "common.h"

int num_procs, for_ser_n;

long long randull(unsigned int *seed)   {
    return ((long long)rand_r(seed) << ((sizeof(int) * 8 - 1) * 2)) | 
           ((long long)rand_r(seed) << ((sizeof(int) * 8 - 1) * 1)) |
           ((long long)rand_r(seed) << ((sizeof(int) * 8 - 1) * 0));
}

void psum(int *data, int n, int *data2) {
  if(n < 2)
    return;

  int i;
  taskfor(i, n/2)
  {
    data2[i] = data[i<<1] + data[(i<<1) + 1];
  } endtaskfor;

  psum(data2, n>>1, data2+(n>>1));

  taskfor(i, n/2)
  {
    data[i<<1] = data2[i]-data[(i<<1)+1];
    data[(i<<1)+1] = data2[i];
  } endtaskfor;

  if(n&1)
    data[n-1] += data2[(n>>1)-1];
}

