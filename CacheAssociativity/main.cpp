#include <iostream>
#include <cstdint>
#include <x86intrin.h>

#define L1d 80*1024
#define L2 1536*1024
#define L3 18*1024*1024

int* createNArray(int numOfFrag, int fragStrideBytes, int cashSizeBytes){
    int fragStrideInt = fragStrideBytes / sizeof(int);

    int arraySize = fragStrideInt * numOfFrag;
    int cashSizeInt = cashSizeBytes/(sizeof (int));
    int arrayFragLen = cashSizeInt/numOfFrag;
    int* array = (int*)calloc(arraySize, sizeof(int));
    if (!array) {
        std::cerr << "allocate error " << arraySize << "\n" ;
        exit(1);
    }


    for(int i = 0; i < arrayFragLen; ++i){
        for(int frag = 0; frag < numOfFrag-1; ++frag){
            array[i + fragStrideInt*frag] = i + fragStrideInt*(frag+1);
        }
        array[i + fragStrideInt * (numOfFrag-1)] = (i+1)%arrayFragLen;
    }

    return array;
}
uint64_t getAverageAccess(const int* array, int N, int round){
    int j=0, l=0;
    for(; j < N; ++j){
        l = array[l];
    }
    uint64_t middle = 0;
    int k=0, i=0;
    uint64_t start = __rdtsc();
    for (; i<N*round; ++i) {
        k = array[k];
    }
    uint64_t end = __rdtsc();

    middle = (end - start)/(N*round);
    return middle;
}

int main() {
    int cacheSize = L3;
    for(int stridePow = 27; stridePow <= 27; ++stridePow){
        int strideBytes = (1 << stridePow);
        std::cout << "strideFrag=" << (1 << stridePow) << "\n";
        for(int numOfFrag = 1; numOfFrag  <= 32; ++numOfFrag ){
            int N;
            auto array = createNArray(numOfFrag , strideBytes, cacheSize);
            auto averageTick = getAverageAccess(array, cacheSize/sizeof(int) , 10);
            std::cout << numOfFrag  << " " << averageTick << "\n";
            free(array);
        }
        std::cout << "\n";
    }
}
