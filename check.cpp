#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

#include "common.h"

int main()  {
    const int n = 1<<17;
    int i;
    bool failed = false;
    dataType *data = new dataType[n];

    for(int k = 0; k < 15; ++k)
    {
        #pragma omp parallel firstprivate(data)
        {
            srand(time(NULL) ^ omp_get_thread_num());
            unsigned int seed = rand();
            #pragma omp for
            for(i = 0; i < n; ++i)  {
                data[i].key = (long long *)randull(&seed);
            }
        }
        failed = false;

        set<long long> ints;
        for(i = 0; i < n; ++i)
            ints.insert((long long)data[i].key);
        
        pSort(data, n, BEST);
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
            return -1;
        }

        printf("Sort succesful.\n");
    }
    return 0;
}
