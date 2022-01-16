#include "../include/rolling_constraint.h"

#include <cmath>

atg_scs::RollingConstraint::RollingConstraint() : Constraint(1, 2) {
    m_local_x = m_local_y = 0.0;
    m_dx = m_dy = 0.0;
    m_radius = 0.0;
    m_ks = 0.0;
    m_kd = 0.0;
}

atg_scs::RollingConstraint::~RollingConstraint() {
    /* void */
}

void atg_scs::RollingConstraint::calculate(
        Output *output,
        int body,
        SystemState *state)
{
    if (body != m_bodies[0]->index && body != m_bodies[1]->index) {
        output->n = 0;
        return;
    }

    output->n = 1;

    const int baseBody = m_bodies[0]->index;
    const int rollingBody = m_bodies[1]->index;

    const double q1_b = state->p_x[baseBody];
    const double q2_b = state->p_y[baseBody];
    const double q3_b = state->theta[baseBody];

    const double q1_r = state->p_x[rollingBody];
    const double q2_r = state->p_y[rollingBody];
    const double q3_r = state->theta[rollingBody];

    const double cos_q3_b = std::cos(q3_b);
    const double sin_q3_b = std::sin(q3_b);

    const double origin_x = q1_b + cos_q3_b * m_local_x - sin_q3_b * m_local_y;
    const double origin_y = q2_b + sin_q3_b * m_local_x + cos_q3_b * m_local_y;
    const double dx = cos_q3_b * m_dx - sin_q3_b * m_dy;
    const double dy = sin_q3_b * m_dx + cos_q3_b * m_dy;

    const double perp_x = -dy;
    const double perp_y = dx;
    
    const double delta_x = q1_r - origin_x;
    const double delta_y = q2_r - origin_y;

    const double s = delta_x * dx + delta_y * dy;

    const double C0 = -q3_r - s * m_radius;
    const double C1 = m_radius - (perp_x * delta_x + perp_y * delta_y);

    if (body == baseBody) {
        const double d_origin_x_dq1 = 1;///
        const double d_origin_x_dq3 = -sin_q3_b * m_local_x - cos_q3_b * m_local_y;///
        
        const double d_origin_y_dq2 = 1;///
        const double d_origin_y_dq3 = cos_q3_b * m_local_x - sin_q3_b * m_local_y;///

        const double d_delta_x_dq1 = -1;///
        const double d_delta_x_dq3 = -d_origin_x_dq3;///

        const double d_delta_y_dq2 = -1;///
        const double d_delta_y_dq3 = -d_origin_y_dq3;///

        const double d_dx_dq3 = -dy;///
        const double d_dy_dq3 = dx;///

        const double ds_dq1 = d_delta_x_dq1 * dx;
        const double ds_dq2 = d_delta_y_dq2 * dy;
        const double ds_dq3 =
            (d_delta_x_dq3 * dx + delta_x * d_dx_dq3) +
            (d_delta_y_dq3 * dy + delta_y * d_dy_dq3);

        // Second derivatives
        const double d2_origin_x_dq3_2 = -d_origin_y_dq3;///
        const double d2_origin_y_dq3_2 = d_origin_x_dq3;///

        const double d2_delta_x_dq3_2 = -d2_origin_x_dq3_2;
        const double d2_delta_y_dq3_2 = -d2_origin_y_dq3_2;

        const double d2_dx_dq3_2 = -d_dy_dq3;
        const double d2_dy_dq3_2 = d_dx_dq3;

        const double d2s_dq1_2 = 0;
        const double d2s_dq2_2 = 0;
        const double d2s_dq3_2 =
            (d2_delta_x_dq3_2 * dx + d_delta_x_dq3 * d_dx_dq3) +
            (d_delta_x_dq3 * d_dx_dq3 + delta_x * d2_dx_dq3_2) +
            (d2_delta_y_dq3_2 * dy + d_delta_y_dq3 * d_dy_dq3) +
            (d_delta_y_dq3 * d_dy_dq3 + delta_y * d2_dy_dq3_2);

        output->dC_dq[0][0] = -ds_dq1 * m_radius;
        output->dC_dq[0][1] = -ds_dq2 * m_radius;
        output->dC_dq[0][2] = -ds_dq3 * m_radius;

        output->d2C_dq2[0][0] = -d2s_dq1_2 * m_radius;
        output->d2C_dq2[0][1] = -d2s_dq2_2 * m_radius;
        output->d2C_dq2[0][2] = -d2s_dq3_2 * m_radius;

        // C1 = m_radius + dy * delta_x - dx * delta_y
        output->dC_dq[0][0] = dy * d_delta_x_dq1;
        output->dC_dq[0][1] = -dx * d_delta_y_dq2;
        output->dC_dq[0][2] =
            (d_dy_dq3 * delta_x + dy * d_delta_x_dq3) -
            (d_dx_dq3 * delta_y + dx * d_delta_y_dq3);

        output->d2C_dq2[0][0] = 0;
        output->d2C_dq2[0][1] = 0;
        output->d2C_dq2[0][2] =
            (d2_dy_dq3_2 * delta_x + d_dy_dq3 * d_delta_x_dq3) +
            (d_dy_dq3 * d_delta_x_dq3 + dy * d2_delta_x_dq3_2) -
            (d2_dx_dq3_2 * delta_y + d_dx_dq3 * d_delta_y_dq3) -
            (d_dx_dq3 * d_delta_y_dq3 + dx * d2_delta_y_dq3_2);
    }
    else if (body == rollingBody) {
        const double d_delta_x_dq1 = 1;
        const double d_delta_x_dq2 = 0;
        const double d_delta_x_dq3 = 0;

        const double d_delta_y_dq1 = 0;
        const double d_delta_y_dq2 = 1;
        const double d_delta_y_dq3 = 0;

        output->dC_dq[0][0] = -1 * dx * m_radius; 
        output->dC_dq[0][1] = -1 * dy * m_radius;
        output->dC_dq[0][2] = -1;

        output->d2C_dq2[0][0] = 0;
        output->d2C_dq2[0][1] = 0;
        output->d2C_dq2[0][2] = 0;

        output->dC_dq[0][0] = dy * d_delta_x_dq1;
        output->dC_dq[0][1] = -dx * d_delta_y_dq2;
        output->dC_dq[0][2] = 0;

        output->d2C_dq2[0][0] = 0;
        output->d2C_dq2[0][1] = 0;
        output->d2C_dq2[0][2] = 0;
    }

    output->ks[0] = m_ks * C0;
    output->kd[0] = m_kd;

    output->ks[0] = m_ks * C1;
    output->kd[0] = m_kd;
}
