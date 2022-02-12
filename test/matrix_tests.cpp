#include <gtest/gtest.h>

#include "utilities.h"

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

TEST(MatrixTests, MatrixLeftScale) {
    const double m0_data[] = {
        -1.0, 2.0,
        3.0, -4.0 };
    const double scale_data[] = {
        1.0,
        2.0 };
    const double scale_matrix_data[] = {
        1.0, 0.0,
        0.0, 2.0 };

    atg_scs::Matrix m0(2, 2);
    atg_scs::Matrix scale_vector(1, 2);
    atg_scs::Matrix scale_matrix(2, 2);
    atg_scs::Matrix resultReference(1, 2), result(1, 2);

    m0.set(m0_data);
    scale_matrix.set(scale_matrix_data);
    scale_vector.set(scale_data);

    scale_matrix.multiply(m0, &resultReference);
    m0.leftScale(scale_vector, &result);

    compareMatrix(result, resultReference);

    m0.destroy();
    scale_vector.destroy();
    scale_matrix.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(MatrixTests, MatrixRightScale) {
    const double m0_data[] = {
        -1.0, 2.0,
        3.0, -4.0 };
    const double scale_data[] = {
        1.0,
        2.0 };
    const double scale_matrix_data[] = {
        1.0, 0.0,
        0.0, 2.0 };

    atg_scs::Matrix m0(2, 2);
    atg_scs::Matrix scale_vector(1, 2);
    atg_scs::Matrix scale_matrix(2, 2);
    atg_scs::Matrix resultReference(1, 2), result(1, 2);

    m0.set(m0_data);
    scale_matrix.set(scale_matrix_data);
    scale_vector.set(scale_data);

    m0.multiply(scale_matrix, &resultReference);
    m0.rightScale(scale_vector, &result);

    compareMatrix(result, resultReference);

    m0.destroy();
    scale_vector.destroy();
    scale_matrix.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(MatrixTests, MatrixTransposeMultiplication) {
    const double m0_data[] = {
        0.0, 1.0,
        -1.0, 0.0 };
    const double v_data[] = {
        1.0,
        2.0 };

    atg_scs::Matrix m0(2, 2);
    atg_scs::Matrix m0_T(2, 2);
    atg_scs::Matrix v(1, 2);
    atg_scs::Matrix resultReference(1, 2);
    atg_scs::Matrix result(1, 2);

    m0.set(m0_data);
    v.set(v_data);

    m0.transpose(&m0_T);
    m0_T.multiply(v, &resultReference);
    m0.transposeMultiply(v, &result);

    compareMatrix(result, resultReference);

    m0.destroy();
    v.destroy();
    result.destroy();
    resultReference.destroy();
    m0_T.destroy();
}
