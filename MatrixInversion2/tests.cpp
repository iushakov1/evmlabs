#include "matrix.h"
#include <gtest/gtest.h>

TEST(Matrix, test1){
    float Mdata[] = {1, 0, 0, 1};
    int Mrow = 2;
    int Mcolumn = 2;
    Matrix M(Mdata, Mrow, Mcolumn);
    auto M_inv = Matrix::getInversed(M, 100);

    int N;
    float expected_data[] = {1, 0, 0, 1};
    float* data = M_inv.getData(N);
    for(int i = 0; i < N; ++i){
        EXPECT_EQ(data[i], expected_data[i]);
    }
    free(data);
}

TEST(Matrix, test2){
    float Mdata[] = {2, 3, 1, 4, 5, 6, 7, 8, 9};
    int Mrow = 3;
    int Mcolumn = 3;
    Matrix M(Mdata, Mrow, Mcolumn);
    auto M_inv = Matrix::getInversed(M, 50000);

    int N;
    float expected_data[] = {(float)(-3)/9, (float)-(19)/9, (float)(13)/9, (float)6/9, (float)(11)/9, (float)(-8)/9, (float)(-3)/9, (float)(5)/9, (float)(-2)/9};
    float* data = M_inv.getData(N);
    for(int i = 0; i < N; ++i){
        if(data[i] > 0){
            EXPECT_TRUE(expected_data[i] - 0.01 <= data[i] && data[i] <= expected_data[i] + 0.01);
        }
        else{
            EXPECT_TRUE((expected_data[i] - 0.01 <= data[i]) && (data[i] <= expected_data[i] + 0.01));
        }
    }
    free(data);
}
