#include "matrix.h"
#include <cstring>
#include <iostream>
#include <xmmintrin.h>

Matrix::Matrix(float *data, int row, int column){
    if(data == nullptr){
        std::cerr << "wrong data ptr for matrix constructor" << std::endl;
        exit(1);
    }
    if(row < 0){
        std::cerr << "wrong row for matrix constructor" << std::endl;
        exit(1);
    }
    if(column < 0){
        std::cerr << "wrong column for matrix constructor" << std::endl;
        exit(1);
    }

    auto* newData = (float*)calloc(row*column, sizeof(float));
    if(newData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }
    std::memcpy(newData, data, row*column*sizeof(float));
    this->data = newData;
    this->row = row;
    this->column = column;
}

Matrix::~Matrix() {
    free(this->data);
}

Matrix Matrix::operator+(Matrix& other) const{
    if (row != other.row || column != other.column) {
        std::cerr << "both matrix must have the same size" << std::endl;
        exit(1);
    }

    int n = row * column;
    auto* resData = (float*)calloc(n, sizeof(float));
    if (!resData) {
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 a = _mm_loadu_ps(&data[i]);
        __m128 b = _mm_loadu_ps(&other.data[i]);
        __m128 r = _mm_add_ps(a, b);
        _mm_storeu_ps(&resData[i], r);
    }

    for (; i < n; ++i) {
        resData[i] = data[i] + other.data[i];
    }

    Matrix res(resData, row, column);
    free(resData);
    return res;
}


Matrix Matrix::operator-(Matrix& other) const {
    if (row != other.row || column != other.column) {
        std::cerr << "both matrix must have the same size" << std::endl;
        exit(1);
    }

    int n = row * column;
    auto* resData = (float*)calloc(n, sizeof(float));
    if (!resData) {
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 a = _mm_loadu_ps(&data[i]);
        __m128 b = _mm_loadu_ps(&other.data[i]);
        __m128 r = _mm_sub_ps(a, b);
        _mm_storeu_ps(&resData[i], r);
    }

    for (; i < n; ++i) {
        resData[i] = data[i] - other.data[i];
    }

    Matrix res(resData, row, column);
    free(resData);
    return res;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (this->column != other.row) {
        std::cerr << "for multiply 1st matrix most have the num of column which is equal to 2nd matrix num of row" << std::endl;
        exit(1);
    }
    int M = this->row;
    int K = this->column;
    int N = other.column;

    Matrix tOther = getTransposed(other);

    auto* resData = (float*)calloc(M * N, sizeof(float));
    if (!resData) {
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            __m128 sum = _mm_setzero_ps();
            int k = 0;

            for (; k + v4fSize <= K; k += v4fSize) {
                __m128 a = _mm_loadu_ps(&this->data[i * K + k]);
                __m128 b = _mm_loadu_ps(&tOther.data[j * K + k]);

                sum = _mm_add_ps(sum, _mm_mul_ps(a, b));
            }

            float tmp[v4fSize];
            _mm_storeu_ps(tmp, sum);
            float val = tmp[0] + tmp[1] + tmp[2] + tmp[3];

            for (; k < K; ++k) {
                val += this->data[i * K + k] * tOther.data[j * K + k];
            }

            resData[i * N + j] = val;
        }
    }

    Matrix res(resData, M, N);
    free(resData);
    return res;
}

float Matrix::getElement(int s, int n) const {
    return data[column * s + n];
}

Matrix Matrix::getInversed(Matrix A, int precision) {
    float A1 = 0;
    float A_inf = 0;
    for(int j = 0; j < A.column; ++j){
        float colSum = 0;
        for(int i = 0; i < A.row; ++i){
            float a_ij = A.getElement(i, j);
            if(a_ij < 0) a_ij = -a_ij;
            colSum += a_ij;
        }
        if(colSum > A1) A1 = colSum;
    }
    for(int i = 0; i < A.row; ++i){
        float rowSum = 0;
        for(int j = 0; j < A.column; ++j){
            float a_ij = A.getElement(i, j);
            if(a_ij < 0) a_ij = -a_ij;
            rowSum += a_ij;
        }
        if(rowSum > A_inf) A_inf = rowSum;
    }
    Matrix B(getTransposed(A));
    B.div(A1*A_inf);
    Matrix I(getIdentityMatrix(A.row));
    auto* buf = (float*)malloc(A.row*A.column*sizeof(float));
    if(buf == nullptr){
        std::cout << "allocating error" << std::endl;
        exit(1);
    }
    Matrix AT = getTransposed(A);
    Matrix BA(B);
    BA.mul(A, buf, AT);
    Matrix R(I);
    R.sub(BA);
    Matrix A_inversed(I);
    int pow = 1;
    Matrix poweredR(R);


    Matrix RT = getTransposed(R);
    while(pow <= precision){
        A_inversed.add(poweredR);
        poweredR.mul(R, buf, RT);
        ++pow;
    }
    Matrix BT = Matrix::getTransposed(B);
    A_inversed.mul(B, buf, BT);
    free(buf);
    return A_inversed;

}

Matrix& Matrix::div(float divisor) {
    for(int i = 0; i < row; ++i){
        for(int j=0; j < column; ++j){
            data[column*i + j] = data[column*i + j] / divisor;
        }
    }
    return *this;
}

Matrix::Matrix(Matrix &other) {
    row = other.row;
    column = other.column;
    auto* newData = (float*)malloc(row*column*sizeof(float));
    if(newData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }
    std::memcpy(newData, other.data, row*column*sizeof(float));
    this->data = newData;
}

Matrix Matrix::getTransposed(const Matrix& A){
    int resRow = A.column;
    int resColumn = A.row;
    auto* resData = (float*)malloc(resRow*resColumn*sizeof(float));
    if(resData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    for(int i = 0; i < A.row; ++i){
        for(int j = 0; j < A.column; ++j){
            resData[j*resColumn + i] = A.data[i*A.column + j];
        }
    }
    Matrix res(resData, resRow, resColumn);
    free(resData);
    return res;
}

Matrix Matrix::getIdentityMatrix(int n) {
    if(n < 1){
        std::cerr << "n must be more than 0" << std::endl;
        exit(1);
    }
    int resRow = n;
    int resColumn = n;
    auto* resData = (float*)calloc(resRow*resColumn, sizeof(float));
    if(resData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }
    for(int i = 0; i < n; ++i){
        resData[i*n + i] = 1;
    }
    Matrix res(resData, resRow, resColumn);
    free(resData);
    return res;
}

int Matrix::operator==(Matrix& other) const{
    if((this->row != other.row) || (this->column != this->row)){
        return 0;
    }
    for(int i = 0; i < this->row; ++i){
        for(int j = 0; j < this->column; ++j){
            if(this->getElement(i, j) != other.getElement(i, j)){
                return 0;
            }
        }
    }

    return 1;
}

int Matrix::operator!=(Matrix& other) const {
    return !(*this == other);
}

Matrix& Matrix::operator=(const Matrix& other) {
    if(&other == this){
        return *this;
    }
    this->row = other.row;
    this->column = other.column;
    free(this->data);
    auto* newData = (float*)calloc(this->row*this->column, sizeof(float));
    if(newData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }
    std::memcpy(newData, other.data, row*column*sizeof(float));
    this->data = newData;
    return *this;
}

void Matrix::printMatrix(Matrix &A) {
    for(int i = 0; i < A.row; ++i){
        for(int j = 0; j < A.column; ++j){
            float a_ij = A.getElement(i, j);
            std::cout << a_ij << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

float *Matrix::getData(int &N) {
    N = row*column;
    auto* res = (float* )calloc(N, sizeof(float));
    if(res == nullptr){
        std::cout << "alocating error" << std::endl;
        exit(1);
    }
    memcpy(res, data, N*sizeof(float));
    return res;
}

void Matrix::add(Matrix& other) {
    if (row != other.row || column != other.column) {
        std::cerr << "both matrix must have the same size" << std::endl;
        exit(1);
    }

    int n = row * column;

    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 a = _mm_loadu_ps(&data[i]);
        __m128 b = _mm_loadu_ps(&other.data[i]);
        __m128 r = _mm_add_ps(a, b);
        _mm_storeu_ps(&this->data[i], r);
    }

    for (; i < n; ++i) {
        this->data[i] = data[i] + other.data[i];
    }

}

void Matrix::sub(Matrix& other) {
    if (row != other.row || column != other.column) {
        std::cerr << "both matrix must have the same size" << std::endl;
        exit(1);
    }

    int n = row * column;

    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 a = _mm_loadu_ps(&data[i]);
        __m128 b = _mm_loadu_ps(&other.data[i]);
        __m128 r = _mm_sub_ps(a, b);
        _mm_storeu_ps(&this->data[i], r);
    }

    for (; i < n; ++i) {
        this->data[i] = data[i] - other.data[i];
    }
}

void Matrix::mul(Matrix& other, float* buf, Matrix& tOther) {
    if (this->column != other.row) {
        std::cerr << "for multiply 1st matrix most have the num of column which is equal to 2nd matrix num of row" << std::endl;
        exit(1);
    }

    int M = this->row;
    int K = this->column;
    int N = other.column;

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            v4f sum = _mm_setzero_ps();
            int k = 0;
            for (; k + v4fSize <= K; k += v4fSize) {
                v4f a = _mm_loadu_ps(&this->data[i * K + k]);
                v4f b = _mm_loadu_ps(&tOther.data[j * K + k]);

                sum = _mm_add_ps(sum, _mm_mul_ps(a, b));
            }

            float tmp[v4fSize];
            _mm_storeu_ps(tmp, sum);
            float val = tmp[0] + tmp[1] + tmp[2] + tmp[3];
            for (; k < K; ++k) {
                val += this->data[i * K + k] * tOther.data[j * K + k];
            }
            buf[i * N + j] = val;
        }
    }

    memcpy(this->data, buf, M * N * sizeof(float));
}

