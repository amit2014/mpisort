#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#include "sort.h"
#include "common.h"

const int ser_n = 1<<9;

void merge(dataType *data, int n1, int n2, dataType *res) {
    // printf("Thread %d doing merge.\n", omp_get_thread_num());
    int i = 0, j = n1, k = 0;

    while(i < n1 && j < n2)
        if((long long)data[i].key < (long long)data[j].key)
            res[k++] = data[i++];
        else
            res[k++] = data[j++];
    
    while(i < n1)
        res[k++] = data[i++];
    while(j < n2)
        res[k++] = data[j++];
}

void mSort_helper(dataType *data, int n, dataType *res)   {
    // printf("Thread %d\n", omp_get_thread_num());
    if(n == 1)  {
        res[0] = data[0];
        return;
    }
    if(n == 0) {
        return;
    }
    if(n < 0)   {
        #ifdef DEBUG
            printf("n < 0 in mSort_helper.\n");
        #endif
        exit(1);
    }
    
    #pragma omp task if(n > ser_n) untied
    mSort_helper(res, n/2, data);
 
    #pragma omp task if(n > ser_n) untied
    mSort_helper(res+n/2, n-n/2, data+n/2);

    #pragma omp taskwait
    merge(data, n/2, n, res);

}

void mSort(dataType *data, int n_total)    {
    MPI_Status status;
    int mpirank, mpisize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    int n = n_total;

    // allocate the extra buffer merge sort needs
    dataType *res = new dataType[n];
    omp_set_num_threads(omp_get_num_procs()*4);

    const int a = mpi_seg_start(0, n, mpirank, mpisize);
    int       b = mpi_seg_end  (0, n, mpirank, mpisize);

    // send data to all nodes
    if(mpirank == 0)    {
        for(int i = 1; i < mpisize; ++i)    {
            int start = mpi_seg_start(0, n, i, mpisize),
                end   = mpi_seg_end  (0, n, i, mpisize);
            if(start < end)
                MPI_Send(data+start, end-start, dataType_MPI(), i,
                        8192+i, MPI_COMM_WORLD);
        }
    }
    else    {
        if(a < b)
            MPI_Recv(data+a, b-a, dataType_MPI(), 0,
                    8192+mpirank, MPI_COMM_WORLD, &status);
    }

    // sort only your segment using omp quick sort
    n = b-a;
    data = data+a;
    res = res+a;
    #pragma omp parallel
    {
        int i;

        #pragma omp for
        for (i = 0; i < n; ++i)
            res[i] = data[i];

        #pragma omp master// implicit nowait
        {
            //printf("%d threads\n", omp_get_num_threads());
            mSort_helper(res, n, data);
        }
    }
    data = data-a;
    res = res-a;
    n = n_total;

    // now we have a sorted segment with each node.
    int factor = 1;
    bool active = true;
    do {
        factor *= 2;
        if(mpirank % factor != 0) {
            // hand over data to parent, and deactivate
            int parent = mpirank - (mpirank % factor);
            MPI_Send(&b, 1, MPI_INT, parent, factor+4096, MPI_COMM_WORLD);
            if(a < b)   {
                MPI_Send(data+a, b-a, dataType_MPI(), parent,
                            factor, MPI_COMM_WORLD);
            }
            active = false;
        }
        else if(mpirank+factor/2 < mpisize) {
            // get data from all child, and merge
            int temp;
            MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE,
                        factor+4096, MPI_COMM_WORLD, &status);
            if(b < temp)    {
                MPI_Recv(data+b, temp-b, dataType_MPI(), MPI_ANY_SOURCE,
                            factor, MPI_COMM_WORLD, &status);
                merge(data+a, b-a, temp-a, res+a);
                b = temp;

                // copy back from res to data
                #pragma omp parallel for 
                for(int i = a; i < b; ++i)
                    data[i] = res[i];
            }
        }
    } while(active && divceil(mpisize, factor) > 1);

    delete [] res;
}
