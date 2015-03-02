#include <stdio.h>
#include <omp.h>
using namespace std;

#include "sort.h"
#include "common.h"

extern int num_procs, for_ser_n;

void mSort(dataType *data, int n);
void qSort(dataType *data, int n);
void rSort(dataType *data, int n);
void bSort(dataType *data, int n);

void pSort(dataType *data, int ndata, SortType sorter)	{
        num_procs = omp_get_num_procs();
        for_ser_n = ((1<<9)*num_procs);

	switch(sorter)	{
		case BEST:  bSort(data, ndata); break;
		case MERGE: mSort(data, ndata); break;
		case QUICK: qSort(data, ndata); break;
                case RADIX: rSort(data, ndata); break;

		default: fprintf(stderr, "Not implemented.");
	}
}
