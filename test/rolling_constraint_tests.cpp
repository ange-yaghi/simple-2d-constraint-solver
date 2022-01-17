#include <gtest/gtest.h>

#include "../include/rolling_constraint.h"

#include <random>
#include <iostream>

void verify(atg_scs::SystemState *state, atg_scs::Constraint *constraint) {
    atg_scs::Constraint::Output o0, o1;
    const int n = constraint->m_bodyCount;
    const int m = constraint->m_constraintCount;

    const double d = 0.001;

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

            const double v0 = q[q_i][j / 3];
            constraint->calculate(&o0, state);
            q[q_i][j / 3] += d;
            constraint->calculate(&o1, state);
            q[q_i][j / 3] = v0;

            const double dC = (o1.C[i] - o0.C[i]) / d;
            EXPECT_NEAR(dC, (o0.dC_dq[i][j] + o1.dC_dq[i][j]) / 2, 1E-4);

            for (int k = 0; k < n * 3; ++k) {
                std::cerr << "d/dq " << k + 1 << "\n";

                const double d2C = (o1.dC_dq[i][k] - o0.dC_dq[i][k]) / d;
                EXPECT_NEAR(d2C, (o0.d2C_dq2[k][i][j] + o1.d2C_dq2[k][i][j]) / 2, 1E-4);
            }
        }
    }
}

TEST(RollingConstraintTests, RollingConstraintTest) {
    atg_scs::RollingConstraint constraint;
    atg_scs::SystemState system;
    atg_scs::Constraint::Output o0, o1;

    system.resize(2);

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

        verify(&system, &constraint);
    }

    system.destroy();
}
