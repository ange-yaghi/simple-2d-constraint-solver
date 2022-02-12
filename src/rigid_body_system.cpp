#include "../include/rigid_body_system.h"

#include <assert.h>
#include <chrono>
#include <cmath>

atg_scs::RigidBodySystem::RigidBodySystem() {
    m_sleSolver = nullptr;
    m_odeSolver = nullptr;

    m_odeSolveMicroseconds = new int[ProfilingSamples];
    m_constraintSolveMicroseconds = new int[ProfilingSamples];
    m_forceEvalMicroseconds = new int[ProfilingSamples];
    m_constraintEvalMicroseconds = new int[ProfilingSamples];
    m_frameIndex = 0;

    for (int i = 0; i < ProfilingSamples; ++i) {
        m_odeSolveMicroseconds[i] = -1;
        m_constraintSolveMicroseconds[i] = -1;
        m_forceEvalMicroseconds[i] = -1;
        m_constraintEvalMicroseconds[i] = -1;
    }
}

atg_scs::RigidBodySystem::~RigidBodySystem() {
    delete[] m_odeSolveMicroseconds;
    delete[] m_constraintSolveMicroseconds;
    delete[] m_forceEvalMicroseconds;
    delete[] m_constraintEvalMicroseconds;
}

void atg_scs::RigidBodySystem::initialize(SleSolver *sleSolver, OdeSolver *odeSolver) {
    m_sleSolver = sleSolver;
    m_odeSolver = odeSolver;
}

void atg_scs::RigidBodySystem::reset() {
    m_rigidBodies.clear();
    m_constraints.clear();
    m_forceGenerators.clear();
}

void atg_scs::RigidBodySystem::addRigidBody(RigidBody *body) {
    m_rigidBodies.push_back(body);
    body->index = (int)m_rigidBodies.size() - 1;
}

void atg_scs::RigidBodySystem::removeRigidBody(RigidBody *body) {
    m_rigidBodies[body->index] = m_rigidBodies.back();
    m_rigidBodies[body->index]->index = body->index;
    m_rigidBodies.resize(m_rigidBodies.size() - 1);
}

atg_scs::RigidBody *atg_scs::RigidBodySystem::getRigidBody(int i) {
    assert(i < m_rigidBodies.size());
    return m_rigidBodies[i];
}

void atg_scs::RigidBodySystem::addConstraint(Constraint *constraint) {
    m_constraints.push_back(constraint);
    constraint->m_index = (int)m_constraints.size() - 1;
}

void atg_scs::RigidBodySystem::removeConstraint(Constraint *constraint) {
    m_constraints[constraint->m_index] = m_constraints.back();
    m_constraints[constraint->m_index]->m_index = constraint->m_index;
    m_constraints.resize(m_constraints.size() - 1);
}

void atg_scs::RigidBodySystem::addForceGenerator(ForceGenerator *forceGenerator) {
    m_forceGenerators.push_back(forceGenerator);
    forceGenerator->m_index = (int)m_forceGenerators.size() - 1;
}

void atg_scs::RigidBodySystem::removeForceGenerator(ForceGenerator *forceGenerator) {
    m_forceGenerators[forceGenerator->m_index] = m_forceGenerators.back();
    m_forceGenerators[forceGenerator->m_index]->m_index = forceGenerator->m_index;
    m_forceGenerators.resize(m_forceGenerators.size() - 1);
}

void atg_scs::RigidBodySystem::process(double dt, int steps) {
    const int n = getRigidBodyCount();

    int
        odeSolveTime = 0,
        constraintSolveTime = 0,
        forceEvalTime = 0,
        constraintEvalTime = 0;

    populateSystemState();
    populateMassMatrices();

    for (int i = 0; i < steps; ++i) {
        m_odeSolver->start(&m_state, dt / steps);

        while (true) {
            const bool done = m_odeSolver->step(&m_state);

            int evalTime = 0, solveTime = 0;

            auto s0 = std::chrono::steady_clock::now();
            processForces();
            auto s1 = std::chrono::steady_clock::now();

            processConstraints(&evalTime, &solveTime);

            auto s2 = std::chrono::steady_clock::now();
            m_odeSolver->solve(&m_state);
            auto s3 = std::chrono::steady_clock::now();

            constraintSolveTime += solveTime;
            constraintEvalTime += evalTime;
            odeSolveTime +=
                std::chrono::duration_cast<std::chrono::microseconds>(s3 - s2).count();
            forceEvalTime +=
                std::chrono::duration_cast<std::chrono::microseconds>(s1 - s0).count();

            if (done) break;
        }

        m_odeSolver->end();
    }

    for (int i = 0; i < n; ++i) {
        m_rigidBodies[i]->v_x = m_state.v_x[i];
        m_rigidBodies[i]->v_y = m_state.v_y[i];

        m_rigidBodies[i]->p_x = m_state.p_x[i];
        m_rigidBodies[i]->p_y = m_state.p_y[i];

        m_rigidBodies[i]->v_theta = m_state.v_theta[i];
        m_rigidBodies[i]->theta = m_state.theta[i];
    }

    m_odeSolveMicroseconds[m_frameIndex] = odeSolveTime;
    m_constraintSolveMicroseconds[m_frameIndex] = constraintSolveTime;
    m_forceEvalMicroseconds[m_frameIndex] = forceEvalTime;
    m_constraintEvalMicroseconds[m_frameIndex] = constraintEvalTime;
    m_frameIndex = (m_frameIndex + 1) % ProfilingSamples;
}

int atg_scs::RigidBodySystem::getFullConstraintCount() const {
    int count = 0;
    for (Constraint *constraint: m_constraints) {
        count += constraint->getConstraintCount();
    }

    return count;
}

float atg_scs::RigidBodySystem::findAverage(int *samples) {
    int accum = 0;
    int count = 0;
    for (int i = 0; i < ProfilingSamples; ++i) {
        if (samples[i] != -1) {
            accum += samples[i];
            ++count;
        }
    }

    if (count == 0) return 0;
    else return (float)accum / count;
}

float atg_scs::RigidBodySystem::getOdeSolveMicroseconds() const {
    return findAverage(m_odeSolveMicroseconds);
}

float atg_scs::RigidBodySystem::getConstraintSolveMicroseconds() const {
    return findAverage(m_constraintSolveMicroseconds);
}

float atg_scs::RigidBodySystem::getConstraintEvalMicroseconds() const {
    return findAverage(m_constraintEvalMicroseconds);
}

float atg_scs::RigidBodySystem::getForceEvalMicroseconds() const {
    return findAverage(m_forceEvalMicroseconds);
}

void atg_scs::RigidBodySystem::populateSystemState() {
    const int n = getRigidBodyCount();

    m_state.resize(n);

    for (int i = 0; i < n; ++i) {
        m_state.a_x[i] = 0;
        m_state.a_y[i] = 0;

        m_state.v_x[i] = m_rigidBodies[i]->v_x;
        m_state.v_y[i] = m_rigidBodies[i]->v_y;

        m_state.p_x[i] = m_rigidBodies[i]->p_x;
        m_state.p_y[i] = m_rigidBodies[i]->p_y;

        m_state.a_theta[i] = 0;
        m_state.v_theta[i] = m_rigidBodies[i]->v_theta;
        m_state.theta[i] = m_rigidBodies[i]->theta;

        m_state.m[i] = m_rigidBodies[i]->m;
    }
}

void atg_scs::RigidBodySystem::populateMassMatrices() {
    const int n = getRigidBodyCount();

    m_iv.M.initialize(3 * n, 3 * n);
    m_iv.M_inv.initialize(3 * n, 3 * n);

    for (int i = 0; i < n; ++i) {
        m_iv.M.set(i * 3 + 0, i * 3 + 0, m_rigidBodies[i]->m);
        m_iv.M.set(i * 3 + 1, i * 3 + 1, m_rigidBodies[i]->m);
        m_iv.M.set(i * 3 + 2, i * 3 + 2, m_rigidBodies[i]->I);

        m_iv.M_inv.set(i * 3 + 0, i * 3 + 0, 1 / m_rigidBodies[i]->m);
        m_iv.M_inv.set(i * 3 + 1, i * 3 + 1, 1 / m_rigidBodies[i]->m);
        m_iv.M_inv.set(i * 3 + 2, i * 3 + 2, 1 / m_rigidBodies[i]->I);
     }
}

void atg_scs::RigidBodySystem::processForces() {
    const int n_f = getForceGeneratorCount();
    const int n = getRigidBodyCount();

    for (int i = 0; i < n; ++i) {
        m_state.f_x[i] = 0.0;
        m_state.f_y[i] = 0.0;
        m_state.t[i] = 0.0;
    }

    for (int i = 0; i < n_f; ++i) {
        m_forceGenerators[i]->apply(&m_state);
    }
}

void atg_scs::RigidBodySystem::processConstraints(
    int *evalTime,
    int *solveTime)
{
    *evalTime = -1;
    *solveTime = -1;

    auto s0 = std::chrono::steady_clock::now();

    const int n = getRigidBodyCount();
    const int m_f = getFullConstraintCount();
    const int m = getConstraintCount();

    m_iv.q_dot.resize(1, n * 3);
    for (int i = 0; i < n; ++i) {
        m_iv.q_dot.set(0, i * 3 + 0, m_state.v_x[i]);
        m_iv.q_dot.set(0, i * 3 + 1, m_state.v_y[i]);
        m_iv.q_dot.set(0, i * 3 + 2, m_state.v_theta[i]);
    }

    m_iv.J.initialize(3 * n, m_f, 0.0);
    m_iv.J_dot.initialize(3 * n, m_f, 0.0);
    m_iv.ks.initialize(1, m_f, 0.0);
    m_iv.kd.initialize(1, m_f, 0.0);
    m_iv.C.initialize(1, m_f, 0.0);

    Constraint::Output constraintOutput;
    int c_i = 0;
    for (int j = 0; j < m; ++j) {
        m_constraints[j]->calculate(&constraintOutput, &m_state);
        const int c_n = m_constraints[j]->getConstraintCount();

        for (int k = 0; k < c_n; ++k) {
            for (int i = 0; i < m_constraints[j]->m_bodyCount * 3; ++i) {
                const int index = m_constraints[j]->m_bodies[i / 3]->index;

                if (index == -1) continue;

                m_iv.J.set(index * 3 + (i % 3), c_i + k,
                        constraintOutput.J[k][i]);

                m_iv.J_dot.set(index * 3 + (i % 3), c_i + k,
                        constraintOutput.J_dot[k][i]);

                m_iv.ks.set(0, c_i + k, constraintOutput.ks[k]);
                m_iv.kd.set(0, c_i + k, constraintOutput.kd[k]);
                m_iv.C.set(0, c_i + k, constraintOutput.C[k]);
            }
        }

        c_i += c_n;
    }

    m_iv.J.multiply(m_iv.q_dot, &m_iv.reg0);
    for (int i = 0; i < m_f; ++i) {
        m_iv.kd.set(0, i, m_iv.kd.get(0, i) * m_iv.reg0.get(0, i));
        m_iv.ks.set(0, i, m_iv.ks.get(0, i) * m_iv.C.get(0, i));
    }

    m_iv.J.transpose(&m_iv.J_T);

    m_iv.F_ext.initialize(1, 3 * n, 0.0);
    for (int i = 0; i < n; ++i) {
        m_iv.F_ext.set(0, i * 3 + 0, m_state.f_x[i]);
        m_iv.F_ext.set(0, i * 3 + 1, m_state.f_y[i]);
        m_iv.F_ext.set(0, i * 3 + 2, m_state.t[i]);
    }

    m_iv.J.multiply(m_iv.M_inv, &m_iv.reg2);
    m_iv.reg2.multiply(m_iv.J_T, &m_iv.left);
    m_iv.reg2.multiply(m_iv.F_ext, &m_iv.reg0);

    m_iv.J_dot.multiply(m_iv.q_dot, &m_iv.reg2);
    m_iv.reg2.negate(&m_iv.reg1);

    m_iv.reg1.subtract(m_iv.reg0, &m_iv.reg2);
    m_iv.reg2.subtract(m_iv.ks, &m_iv.reg0);
    m_iv.reg0.subtract(m_iv.kd, &m_iv.right);

    auto s1 = std::chrono::steady_clock::now();

    const bool solvable =
        m_sleSolver->solve(m_iv.left, m_iv.right, &m_iv.lambda, &m_iv.lambda);
    assert(solvable);

    auto s2 = std::chrono::steady_clock::now();

    m_iv.J_T.multiply(m_iv.lambda, &m_iv.F_C);

    for (int i = 0; i < n; ++i) {
        const double invMass = m_iv.M_inv.get(i * 3 + 0, i * 3 + 0);
        const double invInertia = m_iv.M_inv.get(i * 3 + 2, i * 3 + 2);

        const double F_C_x = (m_f > 0)
            ? m_iv.F_C.get(0, i * 3 + 0)
            : 0;
        const double F_C_y = (m_f > 0)
            ? m_iv.F_C.get(0, i * 3 + 1)
            : 0;
        const double F_C_t = (m_f > 0)
            ? m_iv.F_C.get(0, i * 3 + 2)
            : 0;

        m_state.a_x[i] =
            invMass * (F_C_x + m_iv.F_ext.get(0, i * 3 + 0));
        m_state.a_y[i] =
            invMass * (F_C_y + m_iv.F_ext.get(0, i * 3 + 1));
        m_state.a_theta[i] =
            invInertia * (F_C_t + m_iv.F_ext.get(0, i * 3 + 2));
    }

    auto s3 = std::chrono::steady_clock::now();

    *evalTime =
        std::chrono::duration_cast<std::chrono::microseconds>(s1 - s0 + s3 - s2).count();
    *solveTime =
        std::chrono::duration_cast<std::chrono::microseconds>(s2 - s1).count();
}
