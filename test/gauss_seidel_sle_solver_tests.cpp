#include <gtest/gtest.h>

#include "../include/gaussian_elimination_sle_solver.h"

TEST(GaussianEliminationSleSolverTests, GaussianEliminationSleSolverSanity) {
    atg_scs::GaussianEliminationSleSolver solver;
}

TEST(GaussianEliminationSleSolverTests, GaussSeidelSleSolverBasic) {
    atg_scs::GaussianEliminationSleSolver solver;

    const double L_data[] = {
        500.0, 10.0,
        20.0, -600.0 };
    const double R_data[] = {
        50.0,
        100.0 };

    atg_scs::Matrix solution(1, 2);
    atg_scs::Matrix check(1, 2);
    atg_scs::Matrix L(2, 2);
    atg_scs::Matrix R(1, 2);

    L.set(L_data);
    R.set(R_data);

    const bool solvable = solver.solve(L, R, nullptr, &solution);
    EXPECT_TRUE(solvable);

    L.multiply(solution, &check);

    EXPECT_NEAR(check.get(0, 0), R.get(0, 0), 1e-7);
    EXPECT_NEAR(check.get(0, 1), R.get(0, 1), 1e-7);

    solution.destroy();
    check.destroy();
    L.destroy();
    R.destroy();
}

TEST(GaussianEliminationSleSolverTests, GaussianEliminationSleSolver4x4) {
    atg_scs::GaussianEliminationSleSolver solver;

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
    atg_scs::Matrix R(1, 4);

    L.set(L_data);
    R.set(R_data);

    const bool solvable = solver.solve(L, R, nullptr, &solution);
    EXPECT_TRUE(solvable);

    L.multiply(solution, &check);

    EXPECT_NEAR(check.get(0, 0), R.get(0, 0), 1e-7);
    EXPECT_NEAR(check.get(0, 1), R.get(0, 1), 1e-7);
    EXPECT_NEAR(check.get(0, 2), R.get(0, 2), 1e-7);
    EXPECT_NEAR(check.get(0, 3), R.get(0, 3), 1e-7);

    solution.destroy();
    check.destroy();
    L.destroy();
    R.destroy();
}
