#include <bits/stdc++.h>
#include <omp.h>
#ifdef gprofiler
    #include <google/profiler.h>
#endif
using namespace std;

#include "common.h"

int main()  {
    const int n = 1<<28;
    int i;
    dataType *data = new dataType[n];
    
    #pragma omp parallel firstprivate(data)
    {
        srand(time(NULL) ^ omp_get_thread_num());
        unsigned int seed = rand();
        #pragma omp for
        for(i = 0; i < n; ++i)  {
            data[i].key = (long long *)randull(&seed);
        }
    }

    printf("Beginning sorting..\n");
    
    #ifdef gprofiler
      ProfilerStart("/tmp/profile");
    #endif
    double start, end;
    start = omp_get_wtime();
    pSort(data, n, QUICK);
    end = omp_get_wtime();
    #ifdef gprofiler
      ProfilerStop();
    #endif
    for(i = 0; i < n-1; ++i)
        if((long long)data[i].key > (long long)data[i+1].key) {
            printf("Sort failed!\n");
            printf("Time taken: %.2fs\n", end-start);
            return -1;
        }

    printf("Output looks sorted.\n");
    printf("Time taken: %.2fs\n", end-start);
    
    return 0;
}
