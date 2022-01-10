#include "../include/gravity_force_generator.h"

atg_scs::GravityForceGenerator::GravityForceGenerator() {
    m_g = 9.81;
}

atg_scs::GravityForceGenerator::~GravityForceGenerator() {
    /* void */
}

void atg_scs::GravityForceGenerator::apply(SystemState *state) {
    const int n = state->n;

    for (int i = 0; i < n; ++i) {
        state->f_y[i] += -state->m[i] * m_g;
    }
}
