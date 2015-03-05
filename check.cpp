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

    const int n = 32;
    int i;
    bool failed = false;
    int a = mpi_seg_start(0, n, mpirank, mpisize),
        b = mpi_seg_end  (0, n, mpirank, mpisize);
    dataType *data = new dataType[n];

    for(int k = 0; k < 15; ++k)
    {
        #pragma omp parallel firstprivate(data)
        {
            srand(time(NULL) ^ omp_get_thread_num() ^ mpirank);
            unsigned int seed = rand();
            #pragma omp for
            for(i = a; i < b; ++i)  {
                data[i].key = (long long *)randull(&seed);
            }
        }
        failed = false;
        set<long long> ints;

        if(mpirank == 0)    {
            for(i = 1; i < mpisize; ++i)    {
                int start = mpi_seg_start(0, n, i, mpisize),
                    end   = mpi_seg_end  (0, n, i, mpisize);
                if(start < end)
                    MPI_Recv(data+start, end-start, dataType_MPI(), i,
                            MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            }

            for(i = 0; i < n; ++i)
                ints.insert((long long)data[i].key);
        }
        else if(a < b)   {
            MPI_Send(data+a, b-a, dataType_MPI(), 0, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        pSort(data, n, MERGE);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mpirank == 0)    {
            for(i = 0; i < n-1; ++i)
                if((long long)data[i].key > (long long)data[i+1].key) {
                    printf("(Unordered) ");
                    failed = true;
                    break;
                }

            if(!failed) {
                set<long long> rets;
                for(i = 0; i < n; ++i)  {
                    rets.insert((long long)data[i].key);
                }

                vector<long long> v;
                set_symmetric_difference(
                    ints.begin(), ints.end(),
                    rets.begin(), rets.end(),
                    std::back_inserter(v));
                if(v.size() || ints.size() != rets.size())    {
                    printf("(Numbers changed) ");
                    failed = true;
                }
            }

            if(failed)  {
                printf("Sort failed!\n");
                #ifdef DEBUG
                for(i = 0; i < n; ++i)
                    printf("%lld\n", (long long)data[i].key);
                #endif
                MPI_Finalize();
                return -1;
            }

            printf("Sort succesful.\n");
        }
    }

    MPI_Finalize();
    return 0;
}
