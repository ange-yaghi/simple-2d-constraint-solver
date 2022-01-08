#include "../include/system_state.h"

#include "../include/utilities.h"

#include <assert.h>
#include <cstring>

atg_scs::SystemState::SystemState() {
    a_theta = nullptr;
    v_theta = nullptr;
    theta = nullptr;

    a_x = nullptr;
    a_y = nullptr;
    v_x = nullptr;
    v_y = nullptr;
    p_x = nullptr;
    p_y = nullptr;

    f_x = nullptr;
    f_y = nullptr;
    t = nullptr;

    n = 0;
    dt = 0.0;
}

atg_scs::SystemState::~SystemState() {
    assert(this->n == 0);
}

void atg_scs::SystemState::copy(SystemState *state) {
    resize(state->n);

    if (state->n == 0) {
        return; 
    }

    std::memcpy((void *)a_theta, (void *)state->a_theta, sizeof(double) * this->n);    
    std::memcpy((void *)v_theta, (void *)state->v_theta, sizeof(double) * this->n);    
    std::memcpy((void *)theta, (void *)state->theta, sizeof(double) * this->n);    

    std::memcpy((void *)a_x, (void *)state->a_x, sizeof(double) * this->n);    
    std::memcpy((void *)a_y, (void *)state->a_y, sizeof(double) * this->n);    
    std::memcpy((void *)v_x, (void *)state->v_x, sizeof(double) * this->n);    
    std::memcpy((void *)v_y, (void *)state->v_y, sizeof(double) * this->n);    
    std::memcpy((void *)p_x, (void *)state->p_x, sizeof(double) * this->n);    
    std::memcpy((void *)p_y, (void *)state->p_y, sizeof(double) * this->n);    

    std::memcpy((void *)f_x, (void *)state->f_x, sizeof(double) * this->n);    
    std::memcpy((void *)f_y, (void *)state->f_y, sizeof(double) * this->n);    
    std::memcpy((void *)t, (void *)state->t, sizeof(double) * this->n);    
}

void atg_scs::SystemState::resize(int bodyCount) {
    if (this->n >= bodyCount) {
        return;
    }

    destroy();
    
    this->n = bodyCount;

    a_theta = new double[n];
    v_theta = new double[n];
    theta = new double[n];

    a_x = new double[n];
    a_y = new double[n];
    v_x = new double[n];
    v_y = new double[n];
    p_x = new double[n];
    p_y = new double[n];

    f_x = new double[n];
    f_y = new double[n];
    t = new double[n];
}

void atg_scs::SystemState::destroy() {
    if (this->n > 0) {
        freeArray(a_theta);
        freeArray(v_theta);
        freeArray(theta);

        freeArray(a_x);
        freeArray(a_y);
        freeArray(v_x);
        freeArray(v_y);
        freeArray(p_x);
        freeArray(p_y);

        freeArray(f_x);
        freeArray(f_y);
        freeArray(t);
    }

    n = 0;
}
