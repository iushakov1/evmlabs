#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <x86intrin.h>

#define Round 5

double logBase(double x, double base) {
    if (x <= 0.0) {
        throw std::domain_error("log_base: x must be > 0");
    }
    if (base <= 0.0 || base == 1.0) {
        throw std::domain_error("log_base: base must be > 0 and != 1");
    }
    return std::log(x) / std::log(base);
}

std::string dereferencing(int i){
    switch (i) {
       case 1: {
           return "straight";
       }
       case 2:{
           return "reversed";
       }
       case 3:{
           return "random";
       }
       case 4:{
           return "const offset";
       }
       case 5:{
           return "random block";
       }
       default:{
           return "no such reference";
       }
    }
}

struct average{
    int N;
    int mode;
    uint64_t ticks;
};

void getNScopes(int& nMin, int& nMax){
    nMin = ((4*1024) / sizeof (int));
    nMax = ((32 * 1024 * 1024) / sizeof (int));
}

uint64_t readArray(const int* array, int N){

    for(int t = 0, b = 0; t < N; ++t){
        b = array[b];
    }

    uint64_t middle = 0;
    int k=0, i=0;

    uint64_t start = __rdtsc();
    for (; i<N*Round; ++i)
        k = array[k];
    uint64_t end = __rdtsc();

    middle = (end - start)/(Round*N);
    return middle;
}
int* createNArray(int N, int mode){
    int* array = (int*)calloc(N, sizeof(int));
    if(array == nullptr){
        std::cerr << "memory allocate error" << std::endl;
        exit(1);
    }
    switch (mode) {
        case 1:{
            for(int i = 0; i < N; ++i){
                array[i] = (i + 1)%N;
            }
            break;
        }
        case 2:{
            array[0] = N-1;
            for(int i = 1; i < N; ++i){
                array[i] = i-1;
            }
            break;
        }
        case 3:{
            std::vector<int> perm(N, 0);
            for(int i = 0; i < N; ++i){
                perm[i] = i;
            }
            std::mt19937_64 rng(12345);
            std::shuffle(perm.begin(), perm.end(), rng);
            for (int k = 0; k < N; ++k) {
                int from = perm[k];
                int to   = perm[(k + 1) % N];
                array[from] = to;
            }
            break;
        }
        case 4:{
            int offset = 4000;
            offset += (N%offset == 0);
            for(int i = 0; i < N; ++i){
                array[i] = (i + offset)%N;
            }
            break;
        }
        case 5:{
            int blockSize = 1;
            if (N < blockSize) {
                for(int i = 0; i < N; ++i){
                    array[i] = (i + 1)%N;
                }
                break;
            }

            std::vector<int> perm(N/blockSize, 0);
            for(int i = 0; i < N/blockSize; ++i){
                perm[i] = i;
            }
            std::mt19937_64 rng(12345);
            std::shuffle(perm.begin(), perm.end(), rng);
            for (int k = 0; k < N/blockSize; ++k) {
                int from = perm[k];
                int to   = perm[(k + 1) % (N/blockSize)];
                for(int i = from*blockSize; i < from*blockSize + blockSize; ++i){
                    array[i] = i + 1;
                }
                array[from*blockSize + blockSize - 1] = to*blockSize;
            }

            break;
        }
        default:{
            std::cerr << "no such mode" << std::endl;
            exit(1);
        }
    }
    return array;

}

int main() {
    int NMin, NMax;
    getNScopes(NMin, NMax);
    std::cout << NMax << ' ' << NMin <<'\n' << std::endl;
    std::vector<average> calc;
    for(int mode = 1; mode <= 5; ++mode){
        for(int n = NMin; n <= NMax;){
            int* array = createNArray(n, mode);
            uint64_t ticks = readArray(array, n);
            struct average cur{};
            cur.N = n;
            cur.mode = mode;
            cur.ticks = ticks;
            calc.push_back(cur);
            std::cout << n << std::endl;
            n = int(n*1.3);
            free(array);
        }
    }

    int prevMode = 0;
    for(const auto& res : calc){
        std::string strMode = dereferencing(res.mode);
        if(prevMode != res.mode){
            std::cout << '\n';
            prevMode = res.mode;
        }
        std::cout << "MODE: " << strMode << " N: " << res.N*4/(1024) << "KB" << " TICKS: " << res.ticks << std::endl;

    }


}
