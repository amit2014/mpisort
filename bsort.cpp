#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "sort.h"

const int ser_n = 1<<9;

int round2(int n)   {
    int k = 1;
    while(k < n)
        k <<= 1;
    return k>>1;
}

void bMerge(bool up, dataType *data, int n)   {
    if(n <= 1)
        return;
    int i, m = round2(n);
    dataType tmp;
    for(i = 0; i < n-m; ++i)
    {
        if(((long long)data[i].key > (long long)data[i+m].key) == up) {
            tmp = data[i];
            data[i] = data[i+m];
            data[i+m] = tmp;
        }
    }
    bMerge(up, data, m);
    bMerge(up, data+m, n-m);
}


void bSort_helper(bool up, dataType *data, int n)    {
    if(n <= 1)
        return;

    int m = n>>1;
    #pragma omp task if(m > ser_n) untied
    bSort_helper(!up, data, m);
    #pragma omp task if(n-m > ser_n) untied
    bSort_helper(up, data+m, n-m);
    #pragma omp taskwait
    bMerge(up, data, n);
}

void bSort(dataType *data, int n)   {
    omp_set_num_threads(omp_get_num_procs() * 2);
    #pragma omp parallel firstprivate(data, n)
    {
        #pragma omp master
        {
            // printf("%d threads\n", omp_get_num_threads());
            bSort_helper(true, data, n);
        }
    }
}
