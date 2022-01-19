#include "../include/rolling_constraint.h"

#include <cmath>

atg_scs::RollingConstraint::RollingConstraint() : Constraint(2, 2) {
    m_local_x = m_local_y = 0.0;
    m_dx = m_dy = 0.0;
    m_radius = 0.0;
    m_ks = 10.0;
    m_kd = 1.0;
}

atg_scs::RollingConstraint::~RollingConstraint() {
    /* void */
}

void atg_scs::RollingConstraint::calculate(
        Output *output,
        SystemState *state)
{
    const int baseBody = m_bodies[0]->index;
    const int rollingBody = m_bodies[1]->index;

    const double q1 = state->p_x[baseBody];
    const double q2 = state->p_y[baseBody];
    const double q3 = state->theta[baseBody];

    const double q4 = state->p_x[rollingBody];
    const double q5 = state->p_y[rollingBody];
    const double q6 = state->theta[rollingBody];

    const double cos_q3 = std::cos(q3);
    const double sin_q3 = std::sin(q3);

    const double origin_x = q1 + cos_q3 * m_local_x - sin_q3 * m_local_y;
    const double origin_y = q2 + sin_q3 * m_local_x + cos_q3 * m_local_y;
    const double dx = cos_q3 * m_dx - sin_q3 * m_dy;
    const double dy = sin_q3 * m_dx + cos_q3 * m_dy;

    const double perp_x = -dy;
    const double perp_y = dx;

    const double delta_x = q4 - origin_x;
    const double delta_y = q5 - origin_y;

    const double s = delta_x * dx + delta_y * dy;

    const double C0 = -q6 - s * m_radius;
    const double C1 = m_radius - (perp_x * delta_x + perp_y * delta_y);

    const double d_origin_x_dq3 = -sin_q3 * m_local_x - cos_q3 * m_local_y;
    const double d_origin_y_dq3 = cos_q3 * m_local_x - sin_q3 * m_local_y;

    const double d_delta_x_dq1 = -1;
    const double d_delta_x_dq3 = -d_origin_x_dq3;
    const double d_delta_x_dq4 = 1;

    const double d_delta_y_dq2 = -1;
    const double d_delta_y_dq3 = -d_origin_y_dq3;
    const double d_delta_y_dq5 = 1;

    const double d_dx_dq3 = -dy;
    const double d_dy_dq3 = dx;

    const double d2_dx_dq3_2 = -d_dy_dq3;
    const double d2_dy_dq3_2 = d_dx_dq3;

    const double ds_dq1 = d_delta_x_dq1 * dx;
    const double ds_dq2 = d_delta_y_dq2 * dy;
    const double ds_dq3 =
        (d_delta_x_dq3 * dx + delta_x * d_dx_dq3) +
        (d_delta_y_dq3 * dy + delta_y * d_dy_dq3);

    // Second derivatives
    const double d2_origin_x_dq3_2 = -d_origin_y_dq3;
    const double d2_origin_y_dq3_2 = d_origin_x_dq3;

    const double d2_delta_x_dq3_2 = -d2_origin_x_dq3_2;
    const double d2_delta_y_dq3_2 = -d2_origin_y_dq3_2;

    const double d2s_dq3_2 =
        (d2_delta_x_dq3_2 * dx + d_delta_x_dq3 * d_dx_dq3) +
        (d_delta_x_dq3 * d_dx_dq3 + delta_x * d2_dx_dq3_2) +
        (d2_delta_y_dq3_2 * dy + d_delta_y_dq3 * d_dy_dq3) +
        (d_delta_y_dq3 * d_dy_dq3 + delta_y * d2_dy_dq3_2);

    output->dC_dq[0][0] = -ds_dq1 * m_radius;
    output->dC_dq[0][1] = -ds_dq2 * m_radius;
    output->dC_dq[0][2] = -ds_dq3 * m_radius;

    // C1 = m_radius + dy * delta_x - dx * delta_y
    output->dC_dq[1][0] = dy * d_delta_x_dq1;
    output->dC_dq[1][1] = -dx * d_delta_y_dq2;
    output->dC_dq[1][2] =
        (d_dy_dq3 * delta_x + dy * d_delta_x_dq3) -
        (d_dx_dq3 * delta_y + dx * d_delta_y_dq3);

    output->dC_dq[0][3] = -1 * dx * m_radius;
    output->dC_dq[0][4] = -1 * dy * m_radius;
    output->dC_dq[0][5] = -1;

    output->dC_dq[1][3] = dy * d_delta_x_dq4;
    output->dC_dq[1][4] = -dx * d_delta_y_dq5;
    output->dC_dq[1][5] = 0;

    // Second derivatives

    const double d2C_dq3_2 =
        (d2_dy_dq3_2 * delta_x + d_dy_dq3 * d_delta_x_dq3) +
        (d_dy_dq3 * d_delta_x_dq3 + dy * d2_delta_x_dq3_2) -
        (d2_dx_dq3_2 * delta_y + d_dx_dq3 * d_delta_y_dq3) -
        (d_dx_dq3 * d_delta_y_dq3 + dx * d2_delta_y_dq3_2);

    // d/dq1
    output->d2C_dq2[0][0][0] = 0;
    output->d2C_dq2[0][0][1] = 0;
    output->d2C_dq2[0][0][2] =
        -d_delta_x_dq1 * d_dx_dq3 * m_radius;

    output->d2C_dq2[0][1][0] = 0;
    output->d2C_dq2[0][1][1] = 0;
    output->d2C_dq2[0][1][2] =
        d_dy_dq3 * d_delta_x_dq1;

    output->d2C_dq2[0][0][3] = 0;
    output->d2C_dq2[0][0][4] = 0;
    output->d2C_dq2[0][0][5] = 0;

    output->d2C_dq2[0][1][3] = 0;
    output->d2C_dq2[0][1][4] = 0;
    output->d2C_dq2[0][1][5] = 0;

    // d/dq2
    output->d2C_dq2[1][0][0] = 0;
    output->d2C_dq2[1][0][1] = 0;
    output->d2C_dq2[1][0][2] =
        -d_delta_y_dq2 * d_dy_dq3 * m_radius;

    output->d2C_dq2[1][1][0] = 0;
    output->d2C_dq2[1][1][1] = 0;
    output->d2C_dq2[1][1][2] =
        -d_dx_dq3 * d_delta_y_dq2;

    output->d2C_dq2[1][0][3] = 0;
    output->d2C_dq2[1][0][4] = 0;
    output->d2C_dq2[1][0][5] = 0;

    output->d2C_dq2[1][1][3] = 0;
    output->d2C_dq2[1][1][4] = 0;
    output->d2C_dq2[1][1][5] = 0;

    // d/dq3
    output->d2C_dq2[2][0][0] = -d_delta_x_dq1 * d_dx_dq3 * m_radius;
    output->d2C_dq2[2][0][1] = -d_delta_y_dq2 * d_dy_dq3 * m_radius;
    output->d2C_dq2[2][0][2] = -d2s_dq3_2 * m_radius;

    output->d2C_dq2[2][1][0] = d_dy_dq3 * d_delta_x_dq1;
    output->d2C_dq2[2][1][1] = -d_dx_dq3 * d_delta_y_dq2;
    output->d2C_dq2[2][1][2] = d2C_dq3_2;

    output->d2C_dq2[2][0][3] = -1 * d_dx_dq3 * m_radius;
    output->d2C_dq2[2][0][4] = -1 * d_dy_dq3 * m_radius;
    output->d2C_dq2[2][0][5] = 0;

    output->d2C_dq2[2][1][3] = d_dy_dq3 * d_delta_x_dq4;
    output->d2C_dq2[2][1][4] = -d_dx_dq3 * d_delta_y_dq5;
    output->d2C_dq2[2][1][5] = 0;

    // d/dq4
    output->d2C_dq2[3][0][0] = 0;
    output->d2C_dq2[3][0][1] = 0;
    output->d2C_dq2[3][0][2] = -d_delta_x_dq4 * d_dx_dq3 * m_radius;

    output->d2C_dq2[3][1][0] = 0;
    output->d2C_dq2[3][1][1] = 0;
    output->d2C_dq2[3][1][2] = d_dy_dq3 * d_delta_x_dq4;

    output->d2C_dq2[3][0][3] = 0;
    output->d2C_dq2[3][0][4] = 0;
    output->d2C_dq2[3][0][5] = 0;

    output->d2C_dq2[3][1][3] = 0;
    output->d2C_dq2[3][1][4] = 0;
    output->d2C_dq2[3][1][5] = 0;

    // d/dq5
    output->d2C_dq2[4][0][0] = 0;
    output->d2C_dq2[4][0][1] = 0;
    output->d2C_dq2[4][0][2] = -d_delta_y_dq5 * d_dy_dq3 * m_radius;

    output->d2C_dq2[4][1][0] = 0;
    output->d2C_dq2[4][1][1] = 0;
    output->d2C_dq2[4][1][2] = -d_dx_dq3 * d_delta_y_dq5;

    output->d2C_dq2[4][0][3] = 0;
    output->d2C_dq2[4][0][4] = 0;
    output->d2C_dq2[4][0][5] = 0;

    output->d2C_dq2[4][1][3] = 0;
    output->d2C_dq2[4][1][4] = 0;
    output->d2C_dq2[4][1][5] = 0;

    // d/dq6
    output->d2C_dq2[5][0][0] = 0;
    output->d2C_dq2[5][0][1] = 0;
    output->d2C_dq2[5][0][2] = 0;

    output->d2C_dq2[5][1][0] = 0;
    output->d2C_dq2[5][1][1] = 0;
    output->d2C_dq2[5][1][2] = 0;

    output->d2C_dq2[5][0][3] = 0;
    output->d2C_dq2[5][0][4] = 0;
    output->d2C_dq2[5][0][5] = 0;

    output->d2C_dq2[5][1][3] = 0;
    output->d2C_dq2[5][1][4] = 0;
    output->d2C_dq2[5][1][5] = 0;

    output->ks[0] = m_ks * C0;
    output->kd[0] = m_kd;

    output->ks[1] = m_ks * C1;
    output->kd[1] = m_kd;

    output->C[0] = C0;
    output->C[1] = C1;
}
