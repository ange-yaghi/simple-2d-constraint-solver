#include <gtest/gtest.h>

#include "utilities.h"

#include "../include/gauss_seidel_sle_solver.h"

TEST(GaussSeidelSleSolverTests, GaussSeidelSleSolverSanity) {
    atg_scs::GaussSeidelSleSolver solver;
}

/*
TEST(GaussSeidelSleSolverTests, GaussSeidelSleSolverBasic) {
    atg_scs::GaussSeidelSleSolver solver;

    const double L_data[] = {
        500.0, 10.0,
        20.0, -600.0 };
    const double R_data[] = {
        50.0,
        100.0 };

    atg_scs::Matrix solution(1, 2);
    atg_scs::Matrix check(1, 2);
    atg_scs::Matrix L(2, 2);
    atg_scs::Matrix J(2, 2);
    atg_scs::Matrix R(1, 2);
    atg_scs::Matrix s(1, 2, 1.0);

    J.set(L_data);
    R.set(R_data);

    const bool solvable = solver.solve(J, s, R, nullptr, &solution);
    EXPECT_TRUE(solvable);

    JWJ_t(J, s, &L);
    L.multiply(solution, &check);

    EXPECT_NEAR(check.get(0, 0), R.get(0, 0), 1e-7);
    EXPECT_NEAR(check.get(0, 1), R.get(0, 1), 1e-7);

    solution.destroy();
    check.destroy();
    L.destroy();
    R.destroy();
    s.destroy();
    J.destroy();
}

TEST(GaussSeidelSleSolverTests, GaussSeidelSleSolver4x4) {
    atg_scs::GaussSeidelSleSolver solver;

    const double L_data[] = {
        500.0, 2.0, 3.0, 4.0,
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
    atg_scs::Matrix L(4, 4);
    atg_scs::Matrix J(4, 4);
    atg_scs::Matrix R(1, 4);
    atg_scs::Matrix s(1, 4, 1.0);

    J.set(L_data);
    R.set(R_data);

    const bool solvable = solver.solve(J, s, R, nullptr, &solution);
    EXPECT_TRUE(solvable);

    JWJ_t(J, s, &L);
    L.multiply(solution, &check);

    EXPECT_NEAR(check.get(0, 0), R.get(0, 0), 1e-6);
    EXPECT_NEAR(check.get(0, 1), R.get(0, 1), 1e-6);
    EXPECT_NEAR(check.get(0, 2), R.get(0, 2), 1e-6);
    EXPECT_NEAR(check.get(0, 3), R.get(0, 3), 1e-6);

    solution.destroy();
    check.destroy();
    L.destroy();
    R.destroy();
    J.destroy();
    s.destroy();
}
*/