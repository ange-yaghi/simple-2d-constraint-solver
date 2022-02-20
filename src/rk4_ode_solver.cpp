#include "../include/rk4_ode_solver.h"

atg_scs::Rk4OdeSolver::Rk4OdeSolver() {
    m_stage = m_nextStage = RkStage::Undefined;
}

atg_scs::Rk4OdeSolver::~Rk4OdeSolver() {
    m_initialState.destroy();
    m_accumulator.destroy();
}

void atg_scs::Rk4OdeSolver::start(SystemState *initial, double dt) {
    OdeSolver::start(initial, dt);

    m_initialState.copy(initial);
    m_accumulator.copy(initial);

    m_stage = RkStage::Stage_1;
}

bool atg_scs::Rk4OdeSolver::step(SystemState *state) {
    switch (m_stage) {
        case RkStage::Stage_1:
            state->dt = 0.0;
            break;
        case RkStage::Stage_2:
        case RkStage::Stage_3:
            for (int i = 0; i < state->n; ++i) {
                state->v_theta[i] =
                    m_initialState.v_theta[i] + m_dt * state->a_theta[i] / 2.0;
                state->theta[i] =
                    m_initialState.theta[i] + m_dt * state->v_theta[i] / 2.0;
                state->v_x[i] =
                    m_initialState.v_x[i] + m_dt * state->a_x[i] / 2.0;
                state->v_y[i] =
                    m_initialState.v_y[i] + m_dt * state->a_y[i] / 2.0;
                state->p_x[i] =
                    m_initialState.p_x[i] + m_dt * state->v_x[i] / 2.0;
                state->p_y[i] =
                    m_initialState.p_y[i] + m_dt * state->v_y[i] / 2.0;
            }

            state->dt = m_dt / 2.0;
            break;
        case RkStage::Stage_4:
            for (int i = 0; i < state->n; ++i) {
                state->v_theta[i] =
                    m_initialState.v_theta[i] + m_dt * state->a_theta[i];
                state->theta[i] =
                    m_initialState.theta[i] + m_dt * state->v_theta[i];
                state->v_x[i] =
                    m_initialState.v_x[i] + m_dt * state->a_x[i];
                state->v_y[i] =
                    m_initialState.v_y[i] + m_dt * state->a_y[i];
                state->p_x[i] =
                    m_initialState.p_x[i] + m_dt * state->v_x[i];
                state->p_y[i] =
                    m_initialState.p_y[i] + m_dt * state->v_y[i];
            }

            state->dt = m_dt;
            break;
        default:
            break;
    }

    m_nextStage = getNextStage(m_stage);

    return m_nextStage == RkStage::Complete;
}

void atg_scs::Rk4OdeSolver::solve(SystemState *system) {
    double stageWeight = 0.0;
    switch (m_stage) {
        case RkStage::Stage_1: stageWeight = 1.0; break;
        case RkStage::Stage_2: stageWeight = 2.0; break;
        case RkStage::Stage_3: stageWeight = 2.0; break;
        case RkStage::Stage_4: stageWeight = 1.0; break;
        default: stageWeight = 0.0;
    }

    for (int i = 0; i < system->n; ++i) {
        m_accumulator.v_theta[i] += (m_dt / 6.0) * system->a_theta[i] * stageWeight;
        m_accumulator.theta[i] += (m_dt / 6.0) * system->v_theta[i] * stageWeight;
        m_accumulator.v_x[i] += (m_dt / 6.0) * system->a_x[i] * stageWeight;
        m_accumulator.v_y[i] += (m_dt / 6.0) * system->a_y[i] * stageWeight;
        m_accumulator.p_x[i] += (m_dt / 6.0) * system->v_x[i] * stageWeight;
        m_accumulator.p_y[i] += (m_dt / 6.0) * system->v_y[i] * stageWeight;
    }

    for (int i = 0; i < system->n_c; ++i) {
        m_accumulator.r_x[i] += (m_dt / 6.0) * system->r_x[i] * stageWeight;
        m_accumulator.r_y[i] += (m_dt / 6.0) * system->r_y[i] * stageWeight;
        m_accumulator.r_t[i] += (m_dt / 6.0) * system->r_t[i] * stageWeight;
    }

    if (m_stage == RkStage::Stage_4) {
        for (int i = 0; i < system->n; ++i) {
            system->v_theta[i] = m_accumulator.v_theta[i];
            system->theta[i] = m_accumulator.theta[i];
            system->v_x[i] = m_accumulator.v_x[i];
            system->v_y[i] = m_accumulator.v_y[i];
            system->p_x[i] = m_accumulator.p_x[i];
            system->p_y[i] = m_accumulator.p_y[i];
        }

        for (int i = 0; i < system->n_c; ++i) {
            system->r_x[i] = m_accumulator.r_x[i];
            system->r_y[i] = m_accumulator.r_y[i];
            system->r_t[i] = m_accumulator.r_t[i];
        }
    }

    m_stage = m_nextStage;
}

void atg_scs::Rk4OdeSolver::end() {
    OdeSolver::end();

    m_stage = m_nextStage = RkStage::Undefined;
}

atg_scs::Rk4OdeSolver::RkStage atg_scs::Rk4OdeSolver::getNextStage(RkStage stage) {
    switch (stage) {
        case RkStage::Stage_1: return RkStage::Stage_2;
        case RkStage::Stage_2: return RkStage::Stage_3;
        case RkStage::Stage_3: return RkStage::Stage_4;
        case RkStage::Stage_4: return RkStage::Complete;
        default: return RkStage::Undefined;
    }
}
