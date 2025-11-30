//
// Created by ioanushakov on 24.11.2025.
//

#ifndef MATRIXINVERSION_MATRIX_H
#define MATRIXINVERSION_MATRIX_H

#include <xmmintrin.h>

#define v4fSize 4
typedef __m128 v4f;

class Matrix {
public:
    Matrix(float *data, int row, int column);
    Matrix(Matrix& other);
    ~Matrix();
    Matrix& div(float divisor);

    Matrix operator+(Matrix& other) const;
    Matrix operator-(Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    int operator==(Matrix& other) const ;
    int operator!=(Matrix& other) const;
    Matrix& operator=(const Matrix& other);

    void add(Matrix& other);
    void sub(Matrix& other);
    void mul(Matrix& other, float* buf, Matrix& tOther);

    float getElement(int s, int n) const;
    float* getData(int& N);
    static Matrix getTransposed(const Matrix&);
    static Matrix getInversed(Matrix A, int precision);
    static Matrix getIdentityMatrix(int n);
    static void printMatrix(Matrix& A);

private:
    float* data = nullptr;
    int row = 0;
    int  column = 0;
};


#endif
