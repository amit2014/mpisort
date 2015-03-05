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

    const int n = 1<<17;
    int i;
    bool failed = false;
    dataType *data = new dataType[n];

    for(int k = 0; k < 15; ++k)
    {
        failed = false;
        set<long long> ints;

        if(mpirank == 0)    {
            #pragma omp parallel firstprivate(data)
            {
                srand(time(NULL) ^ omp_get_thread_num() ^ mpirank);
                unsigned int seed = rand();
                #pragma omp for
                for(i = 0; i < n; ++i)  {
                    data[i].key = (long long *)randull(&seed);
                }
            }

            for(i = 0; i < n; ++i)
                ints.insert((long long)data[i].key);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        pSort(data, n, MERGE);

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

            MPI_Bcast(&failed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);

            if(failed)  {
                printf("Sort failed!\n");
                #ifdef DEBUG
                for(i = 0; i < n; ++i)
                    printf("%lld\n", (long long)data[i].key);
                #endif
                MPI_Finalize();
                return 0;
            }

            printf("Sort succesful.\n");
        }
        else    {
            MPI_Bcast(&failed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
            if(failed)  {
                MPI_Finalize();
                return 0;
            }
        }
    }

    MPI_Finalize();
    return 0;
}
