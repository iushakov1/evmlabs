#include "matrix.h"

#include <cstring>
#include <iostream>
#include <cblas.h>

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

Matrix::Matrix(int row, int column) {
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
    this->data = newData;
    this->row = row;
    this->column = column;
}

Matrix::~Matrix() {
    free(this->data);
}

Matrix Matrix::operator+(const Matrix& other) const {
    if((this->row != other.row) || (this->column != other.column)){
        std::cerr << "both matrix most have the same size" << std::endl;
        exit(1);
    }
    int resRow = this->row;
    int resColumn = this->column;
    auto* resData = (float*)calloc(row*column, sizeof(float));
    memcpy(resData, this->data, resColumn*resRow*sizeof(float));
    if(resData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    cblas_saxpy(resColumn*resRow, 1.0f, other.data, 1, resData, 1);
    Matrix res(resData, resRow, resColumn);
    free(resData);
    return res;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if((this->row != other.row) || (this->column != other.column)){
        std::cerr << "both matrix most have the same size" << std::endl;
        exit(1);
    }
    int resRow = this->row;
    int resColumn = this->column;
    auto* resData = (float*)calloc(row*column, sizeof(float));
    memcpy(resData, this->data, resColumn*resRow*sizeof(float));
    if(resData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    cblas_saxpy(resColumn*resRow, -1.0f, other.data, 1, resData, 1);
    Matrix res(resData, resRow, resColumn);
    free(resData);
    return res;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if(this->column != other.row){
        std::cerr << "for multiply 1st matrix most have the num of column which is equal to 2nd matrix num of row" << std::endl;
        exit(1);
    }
    int resRow = this->row;
    int resColumn = other.column;

    auto* resData = (float*)calloc(resRow*resColumn, sizeof(float));
    if(resData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    int M = this->row;
    int K = this->column;
    int N = other.column;

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    M, N, K, 1.0f, this->data, K, other.data, N,
    0.0f, resData, N);
    Matrix res(resData, resRow, resColumn);
    free(resData);
    return res;
}

inline float Matrix::getElement(int s, int n)const {
    return data[column * s + n];
}

int Matrix::setElement(int s, int n, float value) {
    if(s > row || s <= 0){
        std::cerr << "wrong element s index for setting " << std::endl;
        return 1;
    }
    if(n > column || n <= 0){
        std::cerr << "wrong element n index for setting" << std::endl;
        return 1;
    }
    data[column * (s-1) + (n-1)] = value;
    return 0;
}

Matrix Matrix::getInversed(const Matrix& A, int precision) {
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
    B.div(A1);
    B.div(A_inf);
    Matrix I(getIdentityMatrix(A.row));

    Matrix R = (I - (B*A));
    Matrix A_inversed(I);
    int pow = 1;
    Matrix poweredR(R);
    auto* buf = (float*)calloc(poweredR.row*poweredR.column, sizeof(float));
    if(buf == nullptr){
        std::cout << "allocating error" << std::endl;
        exit(1);
    }
    while(pow <= precision){
        A_inversed.add(poweredR);
        poweredR.mul(R, buf);
        ++pow;
    }
    A_inversed = A_inversed*B;
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

    auto* newData = (float*)calloc(row*column, sizeof(float));
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
    auto* resData = (float*)calloc(resRow*resColumn, sizeof(float));
    if(resData == nullptr){
        std::cerr << "error on memory allocation for matrix" << std::endl;
        exit(1);
    }

    for(int i = 0; i < A.row; ++i){
        for(int j = 0; j < A.column; ++j){
            resData[j*resRow + i] = A.data[i*A.column + j];
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

int Matrix::operator==(const Matrix& other) const {
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

int Matrix::operator!=(const Matrix& other) const {
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

    if((this->row != other.row) || (this->column != other.column)){
        std::cerr << "both matrix most have the same size" << std::endl;
        exit(1);
    }
    int resRow = this->row;
    int resColumn = this->column;

    cblas_saxpy(resColumn*resRow, 1.0f, other.data, 1, this->data, 1);
}

void Matrix::sub(Matrix& other) {
    if((this->row != other.row) || (this->column != other.column)){
        std::cerr << "both matrix most have the same size" << std::endl;
        exit(1);
    }
    int resRow = this->row;
    int resColumn = this->column;

    cblas_saxpy(resColumn*resRow, -1.0f, other.data, 1, this->data, 1);
}

void Matrix::mul(Matrix& other, float* buf) {
    if(this->column != other.row){
        std::cerr << "for multiply 1st matrix most have the num of column which is equal to 2nd matrix num of row" << std::endl;
        exit(1);
    }
    int resRow = this->row;
    int resColumn = other.column;

    int M = this->row;
    int K = this->column;
    int N = other.column;

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f, this->data, K, other.data, N,
                0.0f, buf, N);
    memcpy(this->data, buf, resRow * resColumn * sizeof(float));
    memset(buf, 0, resRow * resColumn * sizeof(float));
}