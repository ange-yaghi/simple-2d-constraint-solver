#include <gtest/gtest.h>

#include "utilities.h"

#include "../include/gaussian_elimination_sle_solver.h"

TEST(GaussianEliminationSleSolverTests, GaussianEliminationSleSolverSanity) {
    atg_scs::GaussianEliminationSleSolver solver;
}

TEST(GaussianEliminationSleSolverTests, GaussianEliminationSleSolverBasic) {
    atg_scs::GaussianEliminationSleSolver solver;

    const double J_data[] = {
        500.0, 0.0,
        0.0, -600.0 };
    const double R_data[] = {
        50.0,
        100.0 };

    atg_scs::Matrix solution(1, 2);
    atg_scs::Matrix check(1, 2);
    atg_scs::SparseMatrix J;
    atg_scs::Matrix J_mat(2, 2);
    atg_scs::Matrix R(1, 2);
    atg_scs::Matrix s(1, 2, 1.0);
    atg_scs::Matrix L(2, 2);

    J.initialize(2, 2, 1, 1);
    J.setBlock(0, 0, 0);
    J.setBlock(1, 0, 1);
    J.set(0, 0, 0, 500.0);
    J.set(1, 0, 0, -600.0);

    J_mat.set(J_data);
    R.set(R_data);

    const bool solvable = solver.solve(J, s, R, nullptr, &solution);
    EXPECT_TRUE(solvable);
    
    JWJ_t(J_mat, s, &L);
    L.multiply(solution, &check);

    EXPECT_NEAR(check.get(0, 0), R.get(0, 0), 1e-7);
    EXPECT_NEAR(check.get(0, 1), R.get(0, 1), 1e-7);

    solution.destroy();
    check.destroy();
    J.destroy();
    J_mat.destroy();
    R.destroy();
    s.destroy();
    L.destroy();
}

TEST(GaussianEliminationSleSolverTests, GaussianEliminationSleSolver4x4) {
    atg_scs::GaussianEliminationSleSolver solver;

    const double L_data[] = {
        500.0, 2.0, 0.0, 0.0,
        5.0, -700.0, 45.0, 10.0,
        0.0, 5.0, 200.0, 5.0,
        10.0, 20.0, 30.0, 500.0 };
    const double R_data[] = {
        5.0,
        10.0,
        -1.0,
        20.0 };

    atg_scs::Matrix solution(1, 4);
    atg_scs::Matrix check(1, 4);
    atg_scs::Matrix J(4, 4), JW, J_T;
    atg_scs::SparseMatrix J_sparse, JW_sparse;
    atg_scs::Matrix L(4, 4);
    atg_scs::Matrix R(1, 4);
    atg_scs::Matrix s(1, 4, 1.0);
    atg_scs::Matrix temp;

    J.set(L_data);
    R.set(R_data);

    fullToSparse(J, &J_sparse, 2);
    J_sparse.expand(&temp);
    compareMatrix(temp, J);

    J.rightScale(s, &JW);
    J_sparse.rightScale(s, &JW_sparse);
    JW_sparse.expand(&temp);
    compareMatrix(temp, JW);

    J.transpose(&J_T);
    JW.multiply(J_T, &L);
    JW_sparse.multiplyTranspose(J_sparse, &temp);
    compareMatrix(temp, L);

    const bool solvable = solver.solve(J_sparse, s, R, nullptr, &solution);
    EXPECT_TRUE(solvable);

    L.multiply(solution, &check);

    EXPECT_NEAR(check.get(0, 0), R.get(0, 0), 1e-7);
    EXPECT_NEAR(check.get(0, 1), R.get(0, 1), 1e-7);
    EXPECT_NEAR(check.get(0, 2), R.get(0, 2), 1e-7);
    EXPECT_NEAR(check.get(0, 3), R.get(0, 3), 1e-7);

    solution.destroy();
    check.destroy();
    J.destroy();
    JW.destroy();
    L.destroy();
    R.destroy();
    s.destroy();
    J_sparse.destroy();
    JW_sparse.destroy();
    J_T.destroy();
    temp.destroy();
}
