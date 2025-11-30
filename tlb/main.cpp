#include <iostream>
#include <cstdint>
#include <x86intrin.h>

#define PageSize 4096

int* createTLBarray(int numPages, int pageStrideBytes)
{
    int pageStrideInt = pageStrideBytes / sizeof(int);

    int arraySize = pageStrideInt * numPages;

    int* array = (int*)calloc(arraySize, sizeof(int));
    if (!array) {
        std::cerr << "allocate error\n";
        exit(1);
    }

    for (int p = 0; p < numPages; ++p) {
        int idxHere = p * pageStrideInt;
        int idxNext = ((p + 1) % numPages) * pageStrideInt;
        array[idxHere] = idxNext;
    }

    return array;
}

uint64_t getAverageAccess(const int* array, int numPages, int round) {
    int k = 0;
    for (int i = 0; i < numPages; ++i) {
        k = array[k];
    }

    k = 0;
    uint64_t start = __rdtsc();
    for (int i = 0; i < numPages * round; ++i) {
        k = array[k];
    }
    uint64_t end = __rdtsc();
    uint64_t avg = (end - start) / (numPages * round);

    return avg;
}

int main() {
    for(int stridePow = 0; stridePow <= 6; ++stridePow){
        int strideBytes = PageSize * (1 << stridePow);
        std::cout << "stridePages=" << (1 << stridePow) << "\n";
        for(int numOfPages = 1; numOfPages <= 32; ++numOfPages){
            auto array = createTLBarray(numOfPages, strideBytes);
            auto averageTick = getAverageAccess(array, numOfPages, 200);
            std::cout << numOfPages << " " << averageTick << "\n";
            free(array);
        }
        std::cout << "\n";
    }
}
