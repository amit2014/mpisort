#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

#include "common.h"

int main(int argc, char **argv)  {
    MPI_Init(&argc, &argv);
    int mpirank, mpisize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    const int n = (1<<27);
    int i;
    dataType *data = new dataType[n];
    
    if(mpirank == 0)    {
        #pragma omp parallel firstprivate(data)
        {
            srand(time(NULL) ^ omp_get_thread_num());
            unsigned int seed = rand();
            #pragma omp for
            for(i = 0; i < n; ++i)  {
                data[i].key = (long long *)randull(&seed);
            }
        }
        printf("Generated random numbers..\n");
    }
    
    double start, end;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    pSort(data, n, QUICK);
    end = MPI_Wtime();

    if(mpirank == 0)    {
        for(i = 0; i < n-1; ++i)
            if((long long)data[i].key > (long long)data[i+1].key) {
                printf("Sort failed!\n");
                printf("Time taken: %.2fs\n", end-start);
                MPI_Finalize();
                return -1;
            }

        printf("Output looks sorted.\n");
        printf("Time taken: %.2fs\n", end-start);
    }

    MPI_Finalize();
    return 0;
}
