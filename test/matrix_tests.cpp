#include <gtest/gtest.h>

#include "../include/matrix.h"

TEST(MatrixTests, MatrixInitialization) {
    atg_scs::Matrix matrix(10, 5, 0.0); 

    EXPECT_EQ(matrix.getWidth(), 10);
    EXPECT_EQ(matrix.getHeight(), 5);

    matrix.destroy();
}

TEST(MatrixTests, MatrixMultiplication) {
    const double m0_data[] = {
        0.0, 1.0,
        -1.0, 0.0 };
    const double v_data[] = {
        1.0,
        2.0 };

    atg_scs::Matrix m0(2, 2);
    atg_scs::Matrix v(1, 2);
    atg_scs::Matrix result(1, 2);

    m0.set(m0_data);
    v.set(v_data);

    m0.multiply(v, &result);

    EXPECT_EQ(result.get(0, 0), 2.0);
    EXPECT_EQ(result.get(0, 1), -1.0); 

    m0.destroy();
    v.destroy();
    result.destroy();
}
