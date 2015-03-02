#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "sort.h"

const int ser_n = 1<<14;
unsigned int seed;

void qSort_helper(dataType *data, int n)    {
    // printf ("%d\n", n);
    if(n <= 1)
        return;

    int i, j;
    dataType tmp;
    i = rand_r(&seed) % n;
    tmp = data[i];
    data[i] = data[0];
    data[0] = tmp;

    long long pivot = (long long)data[0].key;
    i = 1; j = n-1;
    while(i <= j)   {
        tmp = data[i];
        data[i] = data[j];
        data[j] = tmp;

        while(i <= j &&
            (long long)data[i].key <= pivot)
                i++;
        while(i <= j &&
            pivot < (long long)data[j].key)
                j--;
        // printf("%d %d\n", i, j);
    }
    if(i-j != 1)    {
        printf("ouch!\n");
        exit(-1);
    }
    tmp = data[0];
    data[0] = data[j];
    data[j] = tmp;

    #pragma omp task if(j > ser_n) untied
    qSort_helper(data, j);

    #pragma omp task if(n-i > ser_n) untied
    qSort_helper(data+i, n-i);
}

void qSort(dataType *data, int n)   {
    #pragma omp parallel firstprivate(data, n) private(seed)
    {
        srand(time(NULL) ^ omp_get_thread_num());
        seed = rand();
        #pragma omp master
        {
            // printf("%d threads\n", omp_get_num_threads());
            qSort_helper(data, n);
        }
    }
}
