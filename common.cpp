#include <stdlib.h>
#include <mpi.h>
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

MPI_Datatype dataType_obj;
bool dataTypeInit = false;

MPI_Datatype dataType_MPI()    {
    if(dataTypeInit)
        return dataType_obj;
    dataTypeInit = true;

    int block_lens[] = {1, LOADSIZE, 1};
    MPI_Aint disps[] = {0, sizeof(long long), sizeof(dataType)};
    MPI_Datatype types[] = {MPI_LONG_LONG, MPI_CHAR, MPI_UB};
    MPI_Type_create_struct(2, block_lens, disps, types, &dataType_obj);
    MPI_Type_commit(&dataType_obj);
    return dataType_obj;
}
