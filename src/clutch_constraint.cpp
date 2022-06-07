#include "../include/clutch_constraint.h"

#include <cmath>
#include <cfloat>

atg_scs::CluchConstraint::CluchConstraint() : Constraint(1, 2) {
    m_ks = 10.0;
    m_kd = 1.0;

    m_maxTorque = 1000.0;
    m_minTorque = -1000.0;
}

atg_scs::CluchConstraint::~CluchConstraint() {
    /* void */
}

void atg_scs::CluchConstraint::calculate(
        Output *output,
        SystemState *state)
{
    const int body = m_bodies[0]->index;
    const int linkedBody = m_bodies[1]->index;

    const double q3 = state->theta[body];
    const double q6 = state->theta[linkedBody];

    const double C = q6 - q3;

    output->J[0][0] = 0.0;
    output->J[0][1] = 0.0;
    output->J[0][2] = -1.0;

    output->J[0][3] = 0.0;
    output->J[0][4] = 0.0;
    output->J[0][5] = 1.0;

    output->J_dot[0][0] = 0;
    output->J_dot[0][1] = 0;
    output->J_dot[0][2] = 0;

    output->J_dot[0][3] = 0;
    output->J_dot[0][4] = 0;
    output->J_dot[0][5] = 0;

    output->kd[0] = m_kd;
    output->ks[0] = m_ks;

    output->C[0] = q6 - q3;

    output->v_bias[0] = 0;
}

void atg_scs::CluchConstraint::limit(
    atg_scs::Matrix *lambda,
    atg_scs::SystemState *state)
{
    const int index = state->indexMap[m_index];
    const double calculatedTorque = lambda->get(0, index);
    //lambda->set(0, index, std::fmin(m_maxTorque, std::fmax(m_minTorque, calculatedTorque)));
}
