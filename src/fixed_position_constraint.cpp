#include "../include/fixed_position_constraint.h"

#include <cmath>

atg_scs::FixedPositionConstraint::FixedPositionConstraint() : Constraint(2, 1) {
    m_local_x = m_local_y = 0.0;
    m_world_x = m_world_y = 0.0;
    m_ks = 500.0;
    m_kd = 20.0;
}

atg_scs::FixedPositionConstraint::~FixedPositionConstraint() {
    /* void */
}

void atg_scs::FixedPositionConstraint::calculate(
        Output *output,
        int body,
        SystemState *state)
{
    if (body != m_bodies[0]->index) {
        output->n = 0;
        return;
    }

    output->n = 2;

    const double q1 = state->p_x[body];
    const double q2 = state->p_y[body];
    const double q3 = state->theta[body];

    const double cos_q3 = std::cos(q3);
    const double sin_q3 = std::sin(q3);

    const double current_x = q1 + cos_q3 * m_local_x - sin_q3 * m_local_y;
    const double current_y = q2 + sin_q3 * m_local_x + cos_q3 * m_local_y;

    const double dx_dq1 = 1.0;
    const double dx_dq2 = 0.0;
    const double dx_dq3 = -sin_q3 * m_local_x - cos_q3 * m_local_y;

    const double dy_dq1 = 0.0;
    const double dy_dq2 = 1.0;
    const double dy_dq3 = cos_q3 * m_local_x - sin_q3 * m_local_y;

    const double d2x_dq1_2 = 0.0;
    const double d2x_dq2_2 = 0.0;
    const double d2x_dq3_2 = -cos_q3 * m_local_x + sin_q3 * m_local_y;

    const double d2y_dq1_2 = 0.0;
    const double d2y_dq2_2 = 0.0;
    const double d2y_dq3_2 = -sin_q3 * m_local_x - cos_q3 * m_local_y;

    const double C1 = current_x - m_world_x;
    const double C2 = current_y - m_world_y;

    output->dC_dq[0][0] = dx_dq1;
    output->dC_dq[0][1] = dx_dq2;
    output->dC_dq[0][2] = dx_dq3;

    output->dC_dq[1][0] = dy_dq1;
    output->dC_dq[1][1] = dy_dq2;
    output->dC_dq[1][2] = dy_dq3;

    output->d2C_dq2[0][0] = d2x_dq1_2;
    output->d2C_dq2[0][1] = d2x_dq2_2;
    output->d2C_dq2[0][2] = d2x_dq3_2;

    output->d2C_dq2[1][0] = d2y_dq1_2;
    output->d2C_dq2[1][1] = d2y_dq2_2;
    output->d2C_dq2[1][2] = d2y_dq3_2;

    output->ks[0] = m_ks * C1;
    output->ks[1] = m_ks * C2;

    output->kd[0] = output->kd[1] = m_kd;
}
