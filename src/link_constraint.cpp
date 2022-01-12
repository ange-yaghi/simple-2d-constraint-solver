#include "../include/link_constraint.h"

#include <cmath>

atg_scs::LinkConstraint::LinkConstraint() : Constraint(2, 2) {
    m_local_x_1 = m_local_y_1 = 0.0;
    m_local_x_2 = m_local_y_2 = 0.0;
    m_ks = 10.0;
    m_kd = 1.0;
}

atg_scs::LinkConstraint::~LinkConstraint() {
    /* void */
}

void atg_scs::LinkConstraint::calculate(
        Output *output,
        int body,
        SystemState *state)
{
    if (body != m_bodies[0]->index && body != m_bodies[1]->index) {
        output->n = 0;
        return;
    }

    double localX, localX_p, localY, localY_p;
    int linkedBody_i;

    if (body == m_bodies[0]->index) {
        localX = m_local_x_1;
        localY = m_local_y_1;
        localX_p = m_local_x_2;
        localY_p = m_local_y_2;
        linkedBody_i = m_bodies[1]->index;
    }
    else {
        localX = m_local_x_2;
        localY = m_local_y_2;
        localX_p = m_local_x_1;
        localY_p = m_local_y_1;
        linkedBody_i = m_bodies[0]->index;
    }

    output->n = 2;

    const double q1 = state->p_x[body];
    const double q2 = state->p_y[body];
    const double q3 = state->theta[body];

    const double q1_p = state->p_x[linkedBody_i];
    const double q2_p = state->p_y[linkedBody_i];
    const double q3_p = state->theta[linkedBody_i];

    const double cos_q3 = std::cos(q3);
    const double sin_q3 = std::sin(q3);

    const double cos_q3_p = std::cos(q3_p);
    const double sin_q3_p = std::sin(q3_p);

    const double bodyX = q1 + cos_q3 * localX - sin_q3 * localY;
    const double bodyY = q2 + sin_q3 * localX + cos_q3 * localY;

    const double linkedBodyX = q1_p + cos_q3_p * localX_p - sin_q3_p * localY_p;
    const double linkedBodyY = q2_p + sin_q3_p * localX_p + cos_q3_p * localY_p;

    const double dbodyX_dq1 = 1.0;
    const double dbodyX_dq2 = 0.0;
    const double dbodyX_dq3 = -sin_q3 * localX - cos_q3 * localY;

    const double dbodyY_dq1 = 0.0;
    const double dbodyY_dq2 = 1.0;
    const double dbodyY_dq3 = cos_q3 * localX - sin_q3 * localY;

    const double d2bodyX_dq1_2 = 0.0;
    const double d2bodyX_dq2_2 = 0.0;
    const double d2bodyX_dq3_2 = -cos_q3 * localX + sin_q3 * localY;

    const double d2bodyY_dq1_2 = 0.0;
    const double d2bodyY_dq2_2 = 0.0;
    const double d2bodyY_dq3_2 = -sin_q3 * localX - cos_q3 * localY;

    const double C1 = bodyX - linkedBodyX;
    const double C2 = bodyY - linkedBodyY;

    if (body == m_bodies[0]->index) {
        output->dC_dq[0][0] = dbodyX_dq1;
        output->dC_dq[0][1] = dbodyX_dq2;
        output->dC_dq[0][2] = dbodyX_dq3;

        output->dC_dq[1][0] = dbodyY_dq1;
        output->dC_dq[1][1] = dbodyY_dq2;
        output->dC_dq[1][2] = dbodyY_dq3;

        output->d2C_dq2[0][0] = d2bodyX_dq1_2;
        output->d2C_dq2[0][1] = d2bodyX_dq2_2;
        output->d2C_dq2[0][2] = d2bodyX_dq3_2;

        output->d2C_dq2[1][0] = d2bodyY_dq1_2;
        output->d2C_dq2[1][1] = d2bodyY_dq2_2;
        output->d2C_dq2[1][2] = d2bodyY_dq3_2;

        output->ks[0] = m_ks * C1;
        output->ks[1] = m_ks * C2;
    }
    else {
        output->dC_dq[0][0] = -dbodyX_dq1;
        output->dC_dq[0][1] = -dbodyX_dq2;
        output->dC_dq[0][2] = -dbodyX_dq3;

        output->dC_dq[1][0] = -dbodyY_dq1;
        output->dC_dq[1][1] = -dbodyY_dq2;
        output->dC_dq[1][2] = -dbodyY_dq3;

        output->d2C_dq2[0][0] = -d2bodyX_dq1_2;
        output->d2C_dq2[0][1] = -d2bodyX_dq2_2;
        output->d2C_dq2[0][2] = -d2bodyX_dq3_2;

        output->d2C_dq2[1][0] = -d2bodyY_dq1_2;
        output->d2C_dq2[1][1] = -d2bodyY_dq2_2;
        output->d2C_dq2[1][2] = -d2bodyY_dq3_2;

        output->ks[0] = m_ks * -C1;
        output->ks[1] = m_ks * -C2;
    }

    output->kd[0] = output->kd[1] = m_kd;
}

void atg_scs::LinkConstraint::setLocalPosition1(double x, double y) {
    m_local_x_1 = x;
    m_local_y_1 = y;
}

void atg_scs::LinkConstraint::setLocalPosition2(double x, double y) {
    m_local_x_2 = x;
    m_local_y_2 = y;
}
