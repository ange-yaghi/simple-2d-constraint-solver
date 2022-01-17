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
        SystemState *state)
{
    const int body = m_bodies[0]->index;
    const int linkedBody = m_bodies[1]->index;

    const double q1 = state->p_x[body];
    const double q2 = state->p_y[body];
    const double q3 = state->theta[body];

    const double q4 = state->p_x[linkedBody];
    const double q5 = state->p_y[linkedBody];
    const double q6 = state->theta[linkedBody];

    const double cos_q3 = std::cos(q3);
    const double sin_q3 = std::sin(q3);

    const double cos_q6 = std::cos(q6);
    const double sin_q6 = std::sin(q6);

    const double bodyX = q1 + cos_q3 * m_local_x_1 - sin_q3 * m_local_y_1;
    const double bodyY = q2 + sin_q3 * m_local_x_1 + cos_q3 * m_local_y_1;

    const double linkedBodyX = q4 + cos_q6 * m_local_x_2 - sin_q6 * m_local_y_2;
    const double linkedBodyY = q5 + sin_q6 * m_local_x_2 + cos_q6 * m_local_y_2;

    // Main body
    const double dbodyX_dq1 = 1.0;
    const double dbodyX_dq2 = 0.0;
    const double dbodyX_dq3 = -sin_q3 * m_local_x_1 - cos_q3 * m_local_y_1;

    const double dbodyY_dq1 = 0.0;
    const double dbodyY_dq2 = 1.0;
    const double dbodyY_dq3 = cos_q3 * m_local_x_1 - sin_q3 * m_local_y_1;

    const double d2bodyX_dq3_2 = -cos_q3 * m_local_x_1 + sin_q3 * m_local_y_1;
    const double d2bodyY_dq3_2 = -sin_q3 * m_local_x_1 - cos_q3 * m_local_y_1;

    // Linked Body
    const double dlinkedBodyX_dq4 = 1.0;
    const double dlinkedBodyX_dq5 = 0.0;
    const double dlinkedBodyX_dq6 = -sin_q6 * m_local_x_2 - cos_q6 * m_local_y_2;

    const double dlinkedBodyY_dq4 = 0.0;
    const double dlinkedBodyY_dq5 = 1.0;
    const double dlinkedBodyY_dq6 = cos_q6 * m_local_x_2 - sin_q6 * m_local_y_2;

    const double d2linkedBodyX_dq6_2 = -cos_q6 * m_local_x_2 + sin_q6 * m_local_y_2;
    const double d2linkedBodyY_dq6_2 = -sin_q6 * m_local_x_2 - cos_q6 * m_local_y_2;

    const double C1 = bodyX - linkedBodyX;
    const double C2 = bodyY - linkedBodyY;

    output->dC_dq[0][0] = dbodyX_dq1;
    output->dC_dq[0][1] = dbodyX_dq2;
    output->dC_dq[0][2] = dbodyX_dq3;

    output->dC_dq[1][0] = dbodyY_dq1;
    output->dC_dq[1][1] = dbodyY_dq2;
    output->dC_dq[1][2] = dbodyY_dq3;

    output->dC_dq[0][3] = -dlinkedBodyX_dq4;
    output->dC_dq[0][4] = -dlinkedBodyX_dq5;
    output->dC_dq[0][5] = -dlinkedBodyX_dq6;

    output->dC_dq[1][3] = -dlinkedBodyY_dq4;
    output->dC_dq[1][4] = -dlinkedBodyY_dq5;
    output->dC_dq[1][5] = -dlinkedBodyY_dq6;

    // d/dq1
    output->d2C_dq2[0][0][0] = 0;
    output->d2C_dq2[0][0][1] = 0;
    output->d2C_dq2[0][0][2] = 0;

    output->d2C_dq2[0][1][0] = 0;
    output->d2C_dq2[0][1][1] = 0;
    output->d2C_dq2[0][1][2] = 0;

    output->d2C_dq2[0][0][3] = 0;
    output->d2C_dq2[0][0][4] = 0;
    output->d2C_dq2[0][0][5] = 0;

    output->d2C_dq2[0][1][3] = 0;
    output->d2C_dq2[0][1][4] = 0;
    output->d2C_dq2[0][1][5] = 0;

    // d/dq2
    output->d2C_dq2[1][0][0] = 0;
    output->d2C_dq2[1][0][1] = 0;
    output->d2C_dq2[1][0][2] = 0;

    output->d2C_dq2[1][1][0] = 0;
    output->d2C_dq2[1][1][1] = 0;
    output->d2C_dq2[1][1][2] = 0;

    output->d2C_dq2[1][0][3] = 0;
    output->d2C_dq2[1][0][4] = 0;
    output->d2C_dq2[1][0][5] = 0;

    output->d2C_dq2[1][1][3] = 0;
    output->d2C_dq2[1][1][4] = 0;
    output->d2C_dq2[1][1][5] = 0;

    // d/dq3
    output->d2C_dq2[2][0][0] = 0;
    output->d2C_dq2[2][0][1] = 0;
    output->d2C_dq2[2][0][2] = d2bodyX_dq3_2;

    output->d2C_dq2[2][1][0] = 0;
    output->d2C_dq2[2][1][1] = 0;
    output->d2C_dq2[2][1][2] = d2bodyY_dq3_2;

    output->d2C_dq2[2][0][3] = 0;
    output->d2C_dq2[2][0][4] = 0;
    output->d2C_dq2[2][0][5] = 0;

    output->d2C_dq2[2][1][3] = 0;
    output->d2C_dq2[2][1][4] = 0;
    output->d2C_dq2[2][1][5] = 0;

    // d/dq4
    output->d2C_dq2[3][0][0] = 0;
    output->d2C_dq2[3][0][1] = 0;
    output->d2C_dq2[3][0][2] = 0;

    output->d2C_dq2[3][1][0] = 0;
    output->d2C_dq2[3][1][1] = 0;
    output->d2C_dq2[3][1][2] = 0;

    output->d2C_dq2[3][0][3] = 0;
    output->d2C_dq2[3][0][4] = 0;
    output->d2C_dq2[3][0][5] = 0;

    output->d2C_dq2[3][1][3] = 0;
    output->d2C_dq2[3][1][4] = 0;
    output->d2C_dq2[3][1][5] = 0;

    // d/dq5
    output->d2C_dq2[4][0][0] = 0;
    output->d2C_dq2[4][0][1] = 0;
    output->d2C_dq2[4][0][2] = 0;

    output->d2C_dq2[4][1][0] = 0;
    output->d2C_dq2[4][1][1] = 0;
    output->d2C_dq2[4][1][2] = 0;

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
    output->d2C_dq2[5][0][5] = -d2linkedBodyX_dq6_2;

    output->d2C_dq2[5][1][3] = 0;
    output->d2C_dq2[5][1][4] = 0;
    output->d2C_dq2[5][1][5] = -d2linkedBodyY_dq6_2;

    output->kd[0] = output->kd[1] = m_kd;
    output->ks[0] = m_ks * C1;
    output->ks[1] = m_ks * C2;
}

void atg_scs::LinkConstraint::setLocalPosition1(double x, double y) {
    m_local_x_1 = x;
    m_local_y_1 = y;
}

void atg_scs::LinkConstraint::setLocalPosition2(double x, double y) {
    m_local_x_2 = x;
    m_local_y_2 = y;
}
