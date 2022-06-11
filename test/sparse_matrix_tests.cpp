#include <gtest/gtest.h>

#include "utilities.h"

#include "../include/matrix.h"
#include "../include/sparse_matrix.h"

TEST(SparseMatrixTests, MatrixTransposeMultiplication) {
    const double m_data[] = {
        0.0, 1.0, 2.0,
        -1.0, 0.0, 5.0,
        4.0, 2.0, 0.0 };

    atg_scs::Matrix m(3, 3);
    atg_scs::Matrix m_T(3, 3);
    atg_scs::SparseMatrix<1> s;
    atg_scs::Matrix resultReference(3, 3);
    atg_scs::Matrix result(3, 3);

    s.initialize(3, 3);
    s.setBlock(0, 0, 1);
    s.setBlock(0, 1, 2);
    s.setBlock(1, 0, 0);
    s.setBlock(1, 1, 2);
    s.setBlock(2, 0, 0);
    s.setBlock(2, 1, 1);
    s.set(0, 0, 0, 1.0);
    s.set(0, 1, 0, 2.0);
    s.set(1, 0, 0, -1.0);
    s.set(1, 1, 0, 5.0);
    s.set(2, 0, 0, 4.0);
    s.set(2, 1, 0, 2.0);

    m.set(m_data);
    m.transpose(&m_T);
    m.multiply(m_T, &resultReference);
    
    s.multiplyTranspose(s, &result);

    compareMatrix(result, resultReference);

    m.destroy();
    m_T.destroy();
    s.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(SparseMatrixTests, MatrixTransposeMultiplicationStride2) {
    const double m_data[] = {
        0.0, 0.0, 1.0, 2.0, 2.0, 3.0,
        -1.0, -2.0, 0.0, 0.0, 5.0, 6.0,
        4.0, 5.0, 2.0, 3.0, 0.0, 0.0 };

    atg_scs::Matrix m(6, 3);
    atg_scs::Matrix m_T;
    atg_scs::SparseMatrix<3> s;
    atg_scs::Matrix resultReference(3, 3);
    atg_scs::Matrix result(3, 3);

    m.set(m_data);
    m.transpose(&m_T);
    fullToSparse(m, &s, 3);

    m.multiply(m_T, &resultReference);
    s.multiplyTranspose(s, &result);

    compareMatrix(result, resultReference);

    m.destroy();
    m_T.destroy();
    s.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(SparseMatrixTests, RightScale) {
    const double m_data[] = {
        0.0, 1.0, 2.0,
        -1.0, 0.0, 5.0,
        4.0, 2.0, 0.0 };
    const double scale_data[] = {
        1.0,
        -1.0,
        3.0 };

    atg_scs::Matrix m(3, 3);
    atg_scs::Matrix scale(1, 3);
    atg_scs::SparseMatrix<1> s;
    atg_scs::SparseMatrix<1> s_scaled;
    atg_scs::Matrix resultReference(3, 3);
    atg_scs::Matrix result(3, 3);

    s.initialize(3, 3);
    s.setBlock(0, 0, 1);
    s.setBlock(0, 1, 2);
    s.setBlock(1, 0, 0);
    s.setBlock(1, 1, 2);
    s.setBlock(2, 0, 0);
    s.setBlock(2, 1, 1);
    s.set(0, 0, 0, 1.0);
    s.set(0, 1, 0, 2.0);
    s.set(1, 0, 0, -1.0);
    s.set(1, 1, 0, 5.0);
    s.set(2, 0, 0, 4.0);
    s.set(2, 1, 0, 2.0);

    m.set(m_data);
    m.rightScale(scale, &resultReference);

    s.rightScale(scale, &s_scaled);
    s_scaled.expand(&result);

    compareMatrix(result, resultReference);

    m.destroy();
    scale.destroy();
    s.destroy();
    s_scaled.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(SparseMatrixTests, LeftScale) {
    const double m_data[] = {
        0.0, 1.0, 2.0,
        -1.0, 0.0, 5.0,
        4.0, 2.0, 0.0 };
    const double scale_data[] = {
        1.0,
        -1.0,
        3.0 };

    atg_scs::Matrix m(3, 3);
    atg_scs::Matrix scale(1, 3);
    atg_scs::SparseMatrix<1> s;
    atg_scs::SparseMatrix<1> s_scaled;
    atg_scs::Matrix resultReference(3, 3);
    atg_scs::Matrix result(3, 3);

    s.initialize(3, 3);
    s.setBlock(0, 0, 1);
    s.setBlock(0, 1, 2);
    s.setBlock(1, 0, 0);
    s.setBlock(1, 1, 2);
    s.setBlock(2, 0, 0);
    s.setBlock(2, 1, 1);
    s.set(0, 0, 0, 1.0);
    s.set(0, 1, 0, 2.0);
    s.set(1, 0, 0, -1.0);
    s.set(1, 1, 0, 5.0);
    s.set(2, 0, 0, 4.0);
    s.set(2, 1, 0, 2.0);

    m.set(m_data);
    m.leftScale(scale, &resultReference);

    s.leftScale(scale, &s_scaled);
    s_scaled.expand(&result);

    compareMatrix(result, resultReference);

    m.destroy();
    scale.destroy();
    s.destroy();
    s_scaled.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(SparseMatrixTests, RightScaleStride2) {
    const double m_data[] = {
        0.0, 0.0, 1.0, 2.0, 2.0, 3.0,
        -1.0, -2.0, 0.0, 0.0, 5.0, 6.0,
        4.0, 5.0, 2.0, 3.0, 0.0, 0.0 };
    const double scale_data[] = {
        1.0,
        -1.0,
        3.0,
        4.0,
        10.0,
        11.0 };

    atg_scs::Matrix m(6, 3);
    atg_scs::Matrix scale(1, 6);
    atg_scs::SparseMatrix<2> s;
    atg_scs::SparseMatrix<2> s_scaled;
    atg_scs::Matrix resultReference(3, 3);
    atg_scs::Matrix result(3, 3);

    s.initialize(6, 3);
    s.setBlock(0, 0, 1);
    s.setBlock(0, 1, 2);
    s.setBlock(1, 0, 0);
    s.setBlock(1, 1, 2);
    s.setBlock(2, 0, 0);
    s.setBlock(2, 1, 1);

    s.set(0, 0, 0, 1.0);
    s.set(0, 0, 1, 2.0);
    s.set(0, 1, 0, 2.0);
    s.set(0, 1, 1, 3.0);

    s.set(1, 0, 0, -1.0);
    s.set(1, 0, 1, -2.0);
    s.set(1, 1, 0, 5.0);
    s.set(1, 1, 1, 6.0);

    s.set(2, 0, 0, 4.0);
    s.set(2, 0, 1, 5.0);
    s.set(2, 1, 0, 2.0);
    s.set(2, 1, 1, 3.0);

    m.set(m_data);
    m.rightScale(scale, &resultReference);

    s.rightScale(scale, &s_scaled);
    s_scaled.expand(&result);

    compareMatrix(result, resultReference);

    m.destroy();
    scale.destroy();
    s.destroy();
    s_scaled.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(SparseMatrixTests, SparseMultiplyingFullMatrix) {
    const double m_data[] = {
        0.0, 0.0, 1.0, 2.0, 2.0, 3.0,
        -1.0, -2.0, 0.0, 0.0, 5.0, 6.0,
        4.0, 5.0, 2.0, 3.0, 0.0, 0.0 };
    const double vectorData[] = {
        1.0,
        -1.0,
        3.0,
        4.0,
        10.0,
        11.0 };

    atg_scs::Matrix m(6, 3);
    atg_scs::Matrix v(1, 6);
    atg_scs::SparseMatrix<3> s;
    atg_scs::Matrix resultReference(1, 3);
    atg_scs::Matrix result(1, 3);

    m.set(m_data);
    v.set(vectorData);

    fullToSparse(m, &s, 3);

    m.multiply(v, &resultReference);
    s.multiply(v, &result);

    compareMatrix(result, resultReference);

    m.destroy();
    v.destroy();
    s.destroy();
    result.destroy();
    resultReference.destroy();
}

TEST(SparseMatrixTests, SparseTransposeMultipleVector) {
    const double m_data[] = {
        0.0, 0.0, 1.0, 2.0, 2.0, 3.0,
        -1.0, -2.0, 0.0, 0.0, 5.0, 6.0,
        4.0, 5.0, 2.0, 3.0, 0.0, 0.0 };
    const double vectorData[] = {
        1.0,
        -1.0,
        3.0 };

    atg_scs::Matrix m(6, 3);
    atg_scs::Matrix v(1, 3);
    atg_scs::SparseMatrix<3> s;
    atg_scs::Matrix resultReference(1, 6);
    atg_scs::Matrix result(1, 6);

    m.set(m_data);
    v.set(vectorData);

    fullToSparse(m, &s, 3);

    m.transposeMultiply(v, &resultReference);
    s.transposeMultiplyVector(v, &result);

    compareMatrix(result, resultReference);

    m.destroy();
    v.destroy();
    s.destroy();
    result.destroy();
    resultReference.destroy();
}
