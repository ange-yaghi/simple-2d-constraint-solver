#include <gtest/gtest.h>

#include "../include/euler_ode_solver.h"

void eulerSolve(int steps, double t, atg_scs::SystemState *state) {
    atg_scs::EulerOdeSolver solver;
    const double dt = t / steps;

    for (int i = 0; i < steps; ++i) {
        solver.start(state, dt);

        while (true) {
            const bool complete = solver.step(state);
            solver.solve(state); 

            if (complete) break;
        }

        solver.end();
    }
}

TEST(EulerOdeSolverTests, EulerOdeSolverSanity) {
    atg_scs::EulerOdeSolver solver;
    atg_scs::SystemState state;

    state.resize(1, 1);

    EXPECT_TRUE(solver.step(&state));

    state.destroy();
}

TEST(EulerOdeSolverTests, EulerOdeSolverIntegrate1Step) {
    atg_scs::EulerOdeSolver solver;
    atg_scs::SystemState state;

    state.resize(1, 1);

    state.a_theta[0] = 10.0;
    state.v_theta[0] = 0.0;
    state.theta[0] = 0.0;

    state.a_x[0] = -2.0;
    state.a_y[0] = -10.0;
    state.v_x[0] = 0.0;
    state.v_y[0] = 0.0;
    state.p_x[0] = 0.0;
    state.p_y[0] = 0.0;

    state.f_x[0] = 0.0;
    state.f_y[0] = 0.0;
    state.t[0] = 0.0;

    eulerSolve(1, 0.1, &state);

    EXPECT_NEAR(state.v_x[0], -0.2, 1E-7);
    EXPECT_NEAR(state.v_y[0], -1.0, 1E-7);
    EXPECT_NEAR(state.p_x[0], -0.01, 1E-1);
    EXPECT_NEAR(state.p_y[0], -0.05, 1E-1);
    EXPECT_NEAR(state.v_theta[0], 1.0, 1E-7);
    EXPECT_NEAR(state.theta[0], 0.0, 1E-7);

    state.destroy();
}

TEST(EulerOdeSolverTests, EulerOdeSolverIntegrate100Steps) {
    atg_scs::EulerOdeSolver solver;
    atg_scs::SystemState state;

    state.resize(1, 1);

    state.a_theta[0] = 10.0;
    state.v_theta[0] = 9.0;
    state.theta[0] = 0.0;

    state.a_x[0] = -2.0;
    state.a_y[0] = -10.0;
    state.v_x[0] = 9.0;
    state.v_y[0] = 5.0;
    state.p_x[0] = 0.0;
    state.p_y[0] = 0.0;

    state.f_x[0] = 0.0;
    state.f_y[0] = 0.0;
    state.t[0] = 0.0;

    eulerSolve(100, 10.0, &state);

    EXPECT_NEAR(state.v_x[0], -11.0, 1E-7);
    EXPECT_NEAR(state.v_y[0], -95.0, 1E-7);
    EXPECT_NEAR(state.p_x[0], -10.0, 1E0);
    EXPECT_NEAR(state.p_y[0], -450.0, 10E0);
    EXPECT_NEAR(state.v_theta[0], 100.0, 10E0);
    EXPECT_NEAR(state.theta[0], 550.0, 50E0);

    state.destroy();
}
