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

    const double q1_dot = state->v_x[baseBody];
    const double q2_dot = state->v_y[baseBody];
    const double q3_dot = state->v_theta[baseBody];

    const double q4_dot = state->v_x[rollingBody];
    const double q5_dot = state->v_y[rollingBody];

    const double cos_q3 = std::cos(q3);
    const double sin_q3 = std::sin(q3);

    const double origin_x = q1 + cos_q3 * m_local_x - sin_q3 * m_local_y;
    const double origin_y = q2 + sin_q3 * m_local_x + cos_q3 * m_local_y;
    const double dx = cos_q3 * m_dx - sin_q3 * m_dy;
    const double dy = sin_q3 * m_dx + cos_q3 * m_dy;

    const double dx_dot = -sin_q3 * q3_dot * m_dx - cos_q3 * q3_dot * m_dy;
    const double dy_dot = cos_q3 * q3_dot * m_dx - sin_q3 * q3_dot * m_dy;

    const double perp_x = -dy;
    const double perp_y = dx;

    const double delta_x = q4 - origin_x;
    const double delta_y = q5 - origin_y;

    const double delta_x_dot =
        q4_dot - (q1_dot - sin_q3 * q3_dot * m_local_x - cos_q3 * q3_dot * m_local_y);
    const double delta_y_dot =
        q5_dot - (q2_dot + cos_q3 * q3_dot * m_local_x - sin_q3 * q3_dot * m_local_y);

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

    const double d_dx_dq3_dot = -cos_q3 * q3_dot * m_dx + sin_q3 * q3_dot * m_dy;
    const double d_dy_dq3_dot = -sin_q3 * q3_dot * m_dx - cos_q3 * q3_dot * m_dy;
    const double d_delta_x_dq3_dot =
        cos_q3 * q3_dot * m_local_x - sin_q3 * q3_dot * m_local_y;
    const double d_delta_y_dq3_dot =
        sin_q3 * q3_dot * m_local_x + cos_q3 * q3_dot * m_local_y;

    const double ds_dq1 = d_delta_x_dq1 * dx;
    const double ds_dq2 = d_delta_y_dq2 * dy;
    const double ds_dq3 =
        (d_delta_x_dq3 * dx + delta_x * d_dx_dq3) +
        (d_delta_y_dq3 * dy + delta_y * d_dy_dq3);

    const double ds_dq1_dot =
        d_delta_x_dq1 * dx_dot;
    const double ds_dq2_dot =
        d_delta_y_dq2 * dy_dot;
    const double ds_dq3_dot =
        (d_delta_x_dq3_dot * dx + d_delta_x_dq3 * dx_dot) +
        (delta_x_dot * d_dx_dq3 + delta_x * d_dx_dq3_dot) +
        (d_delta_y_dq3_dot * dy + d_delta_y_dq3 * dy_dot) +
        (delta_y_dot * d_dy_dq3 + delta_y * d_dy_dq3_dot);

    output->J[0][0] = -ds_dq1 * m_radius;
    output->J[0][1] = -ds_dq2 * m_radius;
    output->J[0][2] = -ds_dq3 * m_radius;

    // C1 = m_radius + dy * delta_x - dx * delta_y
    output->J[1][0] = dy * d_delta_x_dq1;
    output->J[1][1] = -dx * d_delta_y_dq2;
    output->J[1][2] =
        (d_dy_dq3 * delta_x + dy * d_delta_x_dq3) -
        (d_dx_dq3 * delta_y + dx * d_delta_y_dq3);

    output->J[0][3] = -1 * dx * m_radius;
    output->J[0][4] = -1 * dy * m_radius;
    output->J[0][5] = -1;

    output->J[1][3] = dy * d_delta_x_dq4;
    output->J[1][4] = -dx * d_delta_y_dq5;
    output->J[1][5] = 0;

    output->J_dot[0][0] = -ds_dq1_dot * m_radius;
    output->J_dot[0][1] = -ds_dq2_dot * m_radius;
    output->J_dot[0][2] = -ds_dq3_dot * m_radius;

    output->J_dot[1][0] = dy_dot * d_delta_x_dq1;
    output->J_dot[1][1] = -dx_dot * d_delta_y_dq2;
    output->J_dot[1][2] =
        (d_dy_dq3_dot * delta_x + d_dy_dq3 * delta_x_dot) +
        (dy_dot * d_delta_x_dq3 + dy * d_delta_x_dq3_dot) -
        (d_dx_dq3_dot * delta_y + d_dx_dq3 * delta_y_dot) -
        (dx_dot * d_delta_y_dq3 + dx * d_delta_y_dq3_dot);

    output->J_dot[0][3] = -1 * dx_dot * m_radius;
    output->J_dot[0][4] = -1 * dy_dot * m_radius;
    output->J_dot[0][5] = 0;

    output->J_dot[1][3] = dy_dot * d_delta_x_dq4;
    output->J_dot[1][4] = -dx_dot * d_delta_y_dq5;
    output->J_dot[1][5] = 0;

    output->ks[0] = 0;
    output->kd[0] = 0;

    output->ks[1] = m_ks;
    output->kd[1] = m_kd;

    output->C[0] = C0;
    output->C[1] = C1;

    output->v_bias[0] = 0;
    output->v_bias[1] = 0;

    noLimits(output);
}
