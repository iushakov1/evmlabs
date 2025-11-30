
#ifndef MATRIXINVERSION_MATRIX_H
#define MATRIXINVERSION_MATRIX_H

class Matrix {
public:
    Matrix(float *data, int row, int column);
    Matrix(int row, int column);
    Matrix(Matrix& other);
    ~Matrix();
    int setElement(int s, int n, float value);
    Matrix& div(float divisor);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    int operator==(const Matrix& other) const;
    int operator!=(const Matrix& other) const;
    Matrix& operator=(const Matrix& other);

    void add(Matrix& other);
    void sub(Matrix& other);
    void mul(Matrix& other, float* buf);

    float getElement(int s, int n) const;
    float* getData(int& N);
    static Matrix getTransposed(const Matrix&);
    static Matrix getInversed(const Matrix& A, int precision);
    static Matrix getIdentityMatrix(int n);
    static void printMatrix(Matrix& A);
private:
    float* data = nullptr;
    int row = 0;
    int  column = 0;
};


#endif
