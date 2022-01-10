#include "../include/rigid_body.h"

atg_scs::RigidBody::RigidBody() {
    index = 0;
    reset();
}

atg_scs::RigidBody::~RigidBody() {
    /* void */
}

double atg_scs::RigidBody::energy() const {
    const double speed_2 = v_x * v_x + v_y * v_y;
    const double E_k = 0.5 * m * speed_2;
    const double E_r = 0.5 * I * v_theta * v_theta;

    return E_k + E_r;
}

void atg_scs::RigidBody::reset() {
    p_x = p_y = 0.0;
    v_x = v_y = 0.0;

    theta = 0.0;
    v_theta = 0.0;

    m = 0.0;
    I = 0.0;
}
