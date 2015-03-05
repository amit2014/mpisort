#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "sort.h"
#include "common.h"

const int ser_n = 1<<14;
unsigned int seed;

int partition(dataType *data, int n)    {
    // printf ("%d\n", n);
    if(n <= 1)
        return 0;

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

    return j;
}

void qSort_helper(dataType *data, int n)    {
    if(n <= 1)
        return;

    int j = partition(data, n);

    #pragma omp task if(j > ser_n) untied
    qSort_helper(data, j);

    #pragma omp task if(n-j-1 > ser_n) untied
    qSort_helper(data+j+1, n-j-1);
}

void qSort_omp(dataType *data, int n)   {
    if(n <= 1)
        return;

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

void qSort(dataType *data, int n_total)   {
    MPI_Status status;
    MPI_Request nodes_req, data_req;
    int mpirank, mpisize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    srand(time(NULL));
    seed = rand();

    int n = n_total;
    int nodes = mpisize;
    if(mpirank != 0)    {
        // get data to sort from 'someone'
        MPI_Irecv(&nodes, 1, MPI_INT, MPI_ANY_SOURCE, 
                4096+mpirank, MPI_COMM_WORLD, &nodes_req);
        MPI_Irecv(data, n, dataType_MPI(), MPI_ANY_SOURCE,
                mpirank, MPI_COMM_WORLD, &data_req);
        MPI_Wait(&data_req, &status);
        MPI_Get_count(&status, dataType_MPI(), &n);
        MPI_Wait(&nodes_req, &status);

        // first node is guaranteed smallest
        data += 1;
        n--;
    }

    // distribute data 0..n to nodes mpirank..(mpirank+nodes)
    while(n > 1 && nodes > 1) {
        int pivot = partition(data, n);
        int high_nodes = (int)(((double)pivot * (double)nodes) / (double)n);
        if(high_nodes == 0)     high_nodes++;
        if(high_nodes >= nodes)  high_nodes = nodes-1;

        // printf("rank: %d, nodes: %d, sending to %d\n", mpirank, nodes, mpirank+high_nodes);
        nodes -= high_nodes;
        MPI_Send(&nodes, 1, MPI_INT, mpirank+high_nodes,
                    mpirank+high_nodes+4096, MPI_COMM_WORLD);
        MPI_Send(data+pivot, n-pivot, dataType_MPI(), mpirank+high_nodes,
                mpirank+high_nodes, MPI_COMM_WORLD);

        nodes = high_nodes;
        n = pivot;
    }
    // sort left out data
    qSort_omp(data, n);

    if(mpirank != 0)    {
        // remember we got a smallest element as well?
        MPI_Send(data-1, n+1, dataType_MPI(), 0, mpirank, MPI_COMM_WORLD);
    }
    else    {
        for(int i = 1; i < mpisize; ++i)    {
            int recvd = 0;
            MPI_Recv(data+n, n_total, dataType_MPI(), i, i,
                        MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, dataType_MPI(), &recvd);
            n += recvd;
        }
    }
}
