#include "matrix.h"
#include <limits>
#include <iostream>
#include <x86intrin.h>
#include <cstdint>

int main(int argc, char* argv[]) {
    int N; int M;
    if(argc == 3){
        N = std::stoi(argv[1]);
        M = std::stoi(argv[2]);
    }
    else{
        N = 2048;
        M = 10;
    }
    Matrix A = Matrix::getIdentityMatrix(N);
    uint64_t time = 0;
    uint64_t start = __rdtsc();
    A = Matrix::getInversed(A, M);
    uint64_t end = __rdtsc();
    time = (end - start);
    std::cout << "time: " << time << std::endl;
}