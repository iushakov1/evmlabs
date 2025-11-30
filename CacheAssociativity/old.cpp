#include <iostream>
#include <cstdint>
#include <x86intrin.h>

#define L1d 80*1024
#define L2 1536*1024
#define L3 18*1024*1024

/*int* createNArray(int& size, int& N, int numOfFrag, int cashMemB){
    int offset = cashMemB;
    int arrayOffset = offset/(int)sizeof(int);
    int fragLen = cashMemB/numOfFrag;
    int arrayFragLen = (fragLen/(int)sizeof(int));

    int arraySize = arrayOffset * numOfFrag;
    N = arrayFragLen*numOfFrag;
    size = arraySize;

    int* array = (int*)calloc(arraySize, sizeof(int));
    if(array == nullptr){
        std::cerr << "allocate error" << std::endl;
        exit(1);
    }
    for(int i = 0; i < arraySize; ++i){
        array[i] = -1;
    }

    for(int i = 0; i < arrayFragLen; ++i){
        for(int frag = 0; frag < numOfFrag-1; ++frag){
            array[i + arrayOffset*frag] = i + arrayOffset*(frag+1);
        }
        array[i + arrayOffset * (numOfFrag-1)] = (i+1)%arrayFragLen;
    }
    return array;
}*/
int* createNArray(int numOfFrag, unsigned long long fragStrideBytes, unsigned long long cashSizeBytes)
{
    int fragStrideInt = fragStrideBytes / sizeof(int);

    unsigned long long arraySize = fragStrideInt * numOfFrag;
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
    /*for(int numOfFrag = 1; numOfFrag <= 32; ++numOfFrag){
        int N = 0;
        int size = 0;
        auto array = createNArray(size, N, numOfFrag, L3);
        auto averageTick = getAverageAccess(array, N, 100);
        std::cout << numOfFrag << " " << averageTick << std::endl;
        free(array);
    }*/
    int cacheSize = L3;
    for(int stridePow = 27; stridePow <= 27; ++stridePow){
        unsigned long long strideBytes = (1 << stridePow);
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
