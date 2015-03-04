#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

#include "common.h"

int main(int argc, char **argv)  {
    MPI_Init(&argc, &argv);
    MPI_Status status;
    int mpirank, mpisize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    const int n = 1<<25;
    int i;
    int a = mpi_seg_start(0, n, mpirank, mpisize),
        b = mpi_seg_end  (0, n, mpirank, mpisize);
    dataType *data = new dataType[n];
    
    #pragma omp parallel firstprivate(data)
    {
        srand(time(NULL) ^ omp_get_thread_num());
        unsigned int seed = rand();
        #pragma omp for
        for(i = a; i < b; ++i)  {
            data[i].key = (long long *)randull(&seed);
        }
    }

    if(mpirank == 0)    {
        printf("Generated random numbers..\n");
        for(i = 1; i < mpisize; ++i)    {
            int start = mpi_seg_start(0, n, i, mpisize),
                end   = mpi_seg_end  (0, n, i, mpisize);
            if(start < end)
                MPI_Recv(data+start, end-start, dataType_MPI(), i,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }

        printf("Collected data, beginning sorting..\n");
    }
    else if(a < b)   {
        MPI_Send(data+a, b-a, dataType_MPI(), 0, 0, MPI_COMM_WORLD);
    }
    
    double start, end;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    pSort(data, n, MERGE);
    end = MPI_Wtime();

    printf("done.\n");
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
