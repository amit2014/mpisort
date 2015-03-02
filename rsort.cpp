#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "sort.h"
#include "common.h"

const int ser_n = 1<<9;
extern int num_procs, for_ser_n;

void rSort_helper(dataType *data, int n, int shift, dataType *buffer, int *bucket, int *buffer2)    {
    // printf ("%d\n", shift);
    if(n <= 1 || shift < 0)
      return;

    int i, invert = (shift == 63) ? 1 : 0;

    taskfor(i, n)  {
      bucket[i] = (((long long)data[i].key >> shift) &1) ^ invert;
    } endtaskfor;

    psum(bucket, n, buffer2);

    int offset = n-bucket[n-1]-1; // no of zeros - 1
    taskfor(i, n)   {
      if((((long long)data[i].key >> shift) &1) ^ invert)
        buffer[offset + bucket[i]] = data[i];
      else
        buffer[i-bucket[i]] = data[i];
    } endtaskfor;

    taskfor(i, n)
    {
      data[i] = buffer[i];
    } endtaskfor;

    i = ++offset; // i = no of zeros
    #pragma omp task if(i > ser_n) untied
    rSort_helper(data, i, shift-1, buffer, bucket, buffer2);

    #pragma omp task if(n-i > ser_n) untied
    rSort_helper(data+i, n-i, shift-1, buffer+i, bucket+i, buffer2+i);
}

void rSort(dataType *data, int n)   {
    // instead of allocating and deallocating buffers all the time, just allocate once and share
    dataType *newdata = new dataType[n];
    int *bucket = new int[n];
    int *buf2 = new int [n];
    #pragma omp parallel firstprivate(data, n)
    {
        #pragma omp master
        {
            // printf("%d threads\n", omp_get_num_threads());
            rSort_helper(data, n, 63, newdata, bucket, buf2);
        }
    }

    delete [] newdata;
    delete [] bucket;
    delete [] buf2;
}
