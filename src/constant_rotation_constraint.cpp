#include "../include/constant_rotation_constraint.h"

#include <cmath>

atg_scs::ConstantRotationConstraint::ConstantRotationConstraint() : Constraint(1, 1) {
    m_rotationSpeed = 0.0;
    m_ks = 10.0;
    m_kd = 1.0;
}

atg_scs::ConstantRotationConstraint::~ConstantRotationConstraint() {
    /* void */
}

void atg_scs::ConstantRotationConstraint::calculate(
        Output *output,
        SystemState *state)
{
    const int body = m_bodies[0]->index;

    const double q1 = state->p_x[body];
    const double q2 = state->p_y[body];
    const double q3 = state->theta[body];

    const double q3_dot = state->v_theta[body];

    const double C1 = q3_dot - m_rotationSpeed;

    output->J[0][0] = 0.0;
    output->J[0][1] = 0.0;
    output->J[0][2] = 0.0;

    output->J_dot[0][0] = 0;
    output->J_dot[0][1] = 0;
    output->J_dot[0][2] = 0;

    output->ks[0] = m_ks;
    output->ks[1] = m_ks;

    output->kd[0] = output->kd[1] = m_kd;

    output->C[0] = C1;
}
