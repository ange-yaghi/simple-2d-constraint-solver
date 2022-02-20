#include <gtest/gtest.h>

#include "../include/rolling_constraint.h"

#include <random>
#include <iostream>

void verify(atg_scs::SystemState *state, atg_scs::Constraint *constraint) {
    atg_scs::Constraint::Output o0, o1;
    const int n = constraint->m_bodyCount;
    const int m = constraint->getConstraintCount();

    const double d = 0.1;
    const double dt = 0.001;

    for (int i = 0; i < m; ++i) {
        std::cerr << "Constraint " << i << "\n";
        for (int j = 0; j < n * 3; ++j) {
            std::cerr << "q" << j + 1 << "\n";

            const int q_i = j % 3;
            double *q[] = {
                state->p_x,
                state->p_y,
                state->theta
            };

            double *q_dot[] = {
                state->v_x,
                state->v_y,
                state->v_theta
            };

            const double v0 = q[q_i][j / 3];
            q_dot[q_i][j / 3] = d;
            constraint->calculate(&o0, state);
            q[q_i][j / 3] += d * dt;
            constraint->calculate(&o1, state);
            q[q_i][j / 3] = v0;
            q_dot[q_i][j / 3] = 0;

            const double dC = (o1.C[i] - o0.C[i]) / (d * dt);
            EXPECT_NEAR(dC, (o0.J[i][j] + o1.J[i][j]) / 2, 1E-4);

            for (int k = 0; k < n * 3; ++k) {
                for (int l = 0; l < m; ++l) {
                    std::cerr << l << ", " << k << "\n";
                    const double J_dot = (o1.J[l][k] - o0.J[l][k]) / dt;
                    EXPECT_NEAR(J_dot, (o0.J_dot[l][k] + o1.J_dot[l][k]) / 2, 1E-4);
                }
            }
        }
    }
}

TEST(RollingConstraintTests, RollingConstraintTest) {
    atg_scs::RollingConstraint constraint;
    atg_scs::SystemState system;

    system.resize(2, 2);

    atg_scs::RigidBody base, rolling;
    base.index = 0;
    rolling.index = 1;
    constraint.setBaseBody(&base);
    constraint.setRollingBody(&rolling);

    constraint.m_dx = 1.0;
    constraint.m_dy = 0.0;
    constraint.m_local_x = 0.0;
    constraint.m_local_y = 0.1;
    constraint.m_radius = 1.0;

    system.p_x[0] = system.p_x[1] = 0;
    system.p_y[0] = system.p_y[1] = 0;
    system.theta[0] = system.theta[1] = 0;

    std::mt19937 rng;
    rng.seed(0);
    std::uniform_real_distribution<double> realDist;

    const double d = 0.001;
    for (int i = 0; i < 10; ++i) {
        system.p_x[0] = (realDist(rng) - 0.5) * 100;
        system.p_x[1] = (realDist(rng) - 0.5) * 100;
        system.p_y[0] = (realDist(rng) - 0.5) * 100;
        system.p_y[1] = (realDist(rng) - 0.5) * 100;
        system.theta[0] = (realDist(rng) - 0.5) * 100;
        system.theta[1] = (realDist(rng) - 0.5) * 100;
        system.v_x[0] = 0;
        system.v_x[1] = 0;
        system.v_y[0] = 0;
        system.v_y[1] = 0;
        system.v_theta[0] = 0;
        system.v_theta[1] = 0;

        verify(&system, &constraint);
    }

    system.destroy();
}
