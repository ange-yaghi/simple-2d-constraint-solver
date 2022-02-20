#include <gtest/gtest.h>

#include "../include/rk4_ode_solver.h"

void rk4Solve(int steps, double t, atg_scs::SystemState *state) {
    atg_scs::Rk4OdeSolver solver;
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

TEST(Rk4OdeSolverTests, Rk4OdeSolverSanity) {
    atg_scs::Rk4OdeSolver solver;
    atg_scs::SystemState state;

    state.resize(1, 1);

    EXPECT_FALSE(solver.step(&state));

    state.destroy();
}

TEST(Rk4OdeSolverTests, Rk4OdeSolverIntegrate1Step) {
    atg_scs::Rk4OdeSolver solver;
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

    rk4Solve(1, 0.1, &state);

    EXPECT_NEAR(state.v_x[0], -0.2, 1E-7);
    EXPECT_NEAR(state.v_y[0], -1.0, 1E-7);
    EXPECT_NEAR(state.p_x[0], -0.01, 1E-7);
    EXPECT_NEAR(state.p_y[0], -0.05, 1E-7);
    EXPECT_NEAR(state.v_theta[0], 1.0, 1E-7);
    EXPECT_NEAR(state.theta[0], 0.05, 1E-7);

    state.destroy();
}

TEST(Rk4OdeSolverTests, Rk4OdeSolverIntegrate100Steps) {
    atg_scs::Rk4OdeSolver solver;
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

    rk4Solve(100, 10.0, &state);

    EXPECT_NEAR(state.v_x[0], -11.0, 1E-7);
    EXPECT_NEAR(state.v_y[0], -95.0, 1E-7);
    EXPECT_NEAR(state.p_x[0], -10.0, 1E-7);
    EXPECT_NEAR(state.p_y[0], -450.0, 1E-7);
    EXPECT_NEAR(state.v_theta[0], 109.0, 1E-7);
    EXPECT_NEAR(state.theta[0], 590.0, 1E-7);

    state.destroy();
}
