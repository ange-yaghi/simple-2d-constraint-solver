#include "../include/fixed_rotation_constraint.h"

atg_scs::FixedRotationConstraint::FixedRotationConstraint() : Constraint(1, 1) {
    m_rotation = 0;
    m_ks = 10.0;
    m_kd = 1.0;
}

atg_scs::FixedRotationConstraint::~FixedRotationConstraint() {
    /* void */
}

void atg_scs::FixedRotationConstraint::calculate(Output *output, SystemState *state) {
    const int body = m_bodies[0]->index;

    const double q3 = state->theta[body];
    const double C = q3 - m_rotation;

    output->J[0][0] = 0;
    output->J[0][1] = 0;
    output->J[0][2] = 1.0;

    output->J_dot[0][0] = 0;
    output->J_dot[0][1] = 0;
    output->J_dot[0][2] = 0;

    output->ks[0] = m_ks;
    output->kd[0] = m_kd;

    output->C[0] = C;

    output->v_bias[0] = 0;

    noLimits(output);
}
