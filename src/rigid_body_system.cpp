#include "../include/rigid_body_system.h"

#include <assert.h>

atg_scs::RigidBodySystem::RigidBodySystem() {
    m_sleSolver = nullptr;
    m_odeSolver = nullptr;
}

atg_scs::RigidBodySystem::~RigidBodySystem() {
    /* void */
}

void atg_scs::RigidBodySystem::initialize(SleSolver *sleSolver, OdeSolver *odeSolver) {
    m_sleSolver = sleSolver;
    m_odeSolver = odeSolver;
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

void atg_scs::RigidBodySystem::process(double dt) {
    const int n = getRigidBodyCount();
    const int m_f = getFullConstraintCount();
    const int m = getConstraintCount();

    populateSystemState();
    populateMassMatrices();

    const int steps = 1;
    for (int i = 0; i < steps; ++i) {
        m_odeSolver->start(&m_state, dt / steps);

        while (true) {
            const bool done = m_odeSolver->step(&m_state);

            processForces();
            processConstraints();

            m_odeSolver->solve(&m_state);

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

    for (int i = 0, c_i = 0; i < m; ++i) {
        const int n_sub = m_constraints[i]->m_constraintCount;
        for (int j = 0; j < n_sub; ++j, ++c_i) {
            for (int k = 0; k < m_constraints[i]->m_bodyCount; ++k) {
                m_constraints[i]->m_f_x[k][j] =
                    m_iv.R.get(c_i, m_constraints[i]->m_bodies[k]->index * 3 + 0);
                m_constraints[i]->m_f_x[k][j] =
                    m_iv.R.get(c_i, m_constraints[i]->m_bodies[k]->index * 3 + 1);
                m_constraints[i]->m_t[k][j] =
                    m_iv.R.get(c_i, m_constraints[i]->m_bodies[k]->index * 3 + 2);
            }
        }
    }
}

int atg_scs::RigidBodySystem::getFullConstraintCount() const {
    int count = 0;
    for (Constraint *constraint: m_constraints) {
        count += constraint->m_constraintCount;
    }

    return count;
}

void atg_scs::RigidBodySystem::populateSystemState() {
    const int n = getRigidBodyCount();

    m_state.resize(n);

    for (int i = 0; i < n; ++i) {
        m_state.v_x[i] = m_rigidBodies[i]->v_x;
        m_state.v_y[i] = m_rigidBodies[i]->v_y;

        m_state.p_x[i] = m_rigidBodies[i]->p_x;
        m_state.p_y[i] = m_rigidBodies[i]->p_y;

        m_state.v_theta[i] = m_rigidBodies[i]->v_theta;
        m_state.theta[i] = m_rigidBodies[i]->theta;
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

void atg_scs::RigidBodySystem::processConstraints() {
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
    m_iv.C_ks.initialize(1, m_f, 0.0);
    m_iv.C_kd.initialize(1, m_f, 0.0);

    Constraint::Output constraintOutput;
    int c_i = 0;
    for (int j = 0;  j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            m_constraints[j]->calculate(&constraintOutput, i, &m_state);

            for (int k = 0; k < constraintOutput.n; ++k) {
                m_iv.J.set(i * 3 + 0, c_i + k, constraintOutput.dC_dq[k][0]);
                m_iv.J.set(i * 3 + 1, c_i + k, constraintOutput.dC_dq[k][1]);
                m_iv.J.set(i * 3 + 2, c_i + k, constraintOutput.dC_dq[k][2]);

                m_iv.J_dot.set(i * 3 + 0, c_i + k,
                    constraintOutput.d2C_dq2[k][0] * m_state.v_x[i]);
                m_iv.J_dot.set(i * 3 + 1, c_i + k,
                    constraintOutput.d2C_dq2[k][1] * m_state.v_y[i]);
                m_iv.J_dot.set(i * 3 + 2, c_i + k,
                    constraintOutput.d2C_dq2[k][2] * m_state.v_theta[i]);

                m_iv.C_ks.set(0, c_i + k, constraintOutput.ks[k]);
                m_iv.C_kd.set(0, c_i + k, constraintOutput.kd[k]);
            }
        }

        c_i += m_constraints[j]->m_constraintCount;
    }

    m_iv.J.multiply(m_iv.q_dot, &m_iv.reg0);
    for (int i = 0; i < m_f; ++i) {
        m_iv.C_kd.set(0, i, m_iv.C_kd.get(0, i) * m_iv.reg0.get(0, i));
    }

    m_iv.J.transpose(&m_iv.J_T);

    m_iv.F_ext.initialize(1, 3 * n, 0.0);
    for (int i = 0; i < n; ++i) {
        m_iv.F_ext.set(0, i * 3 + 0, m_state.f_x[i]);
        m_iv.F_ext.set(0, i * 3 + 1, m_state.f_y[i]);
        m_iv.F_ext.set(0, i * 3 + 2, m_state.t[i]);
    }

    m_iv.J.multiply(m_iv.M_inv, &m_iv.reg0);
    m_iv.reg0.multiply(m_iv.J_T, &m_iv.left);

    m_iv.J_dot.multiply(m_iv.q_dot, &m_iv.reg0);
    m_iv.reg0.negate(&m_iv.reg1);

    m_iv.J.multiply(m_iv.M_inv, &m_iv.reg2);
    m_iv.reg2.multiply(m_iv.F_ext, &m_iv.reg0);

    m_iv.reg1.subtract(m_iv.reg0, &m_iv.reg2);
    m_iv.reg2.subtract(m_iv.C_ks, &m_iv.reg1);
    m_iv.reg1.subtract(m_iv.C_kd, &m_iv.right);

    const bool solvable =
        m_sleSolver->solve(m_iv.left, m_iv.right, &m_iv.lambda, &m_iv.lambda);
    assert(solvable);

    m_iv.J_T.multiply(m_iv.lambda, &m_iv.F_C);

    m_iv.lambdaScale.initialize(m_f, m_f, 0.0);
    for (int i = 0; i < m_f; ++i) {
        m_iv.lambdaScale.set(i, i, m_iv.lambda.get(0, i));
    }

    m_iv.J_T.multiply(m_iv.lambdaScale, &m_iv.R);

    for (int i = 0; i < n; ++i) {
        const double invMass = m_iv.M_inv.get(i * 3 + 0, i * 3 + 0);
        const double invInertia = m_iv.M_inv.get(i * 3 + 2, i * 3 + 2);

        m_state.a_x[i] =
            invMass * (m_iv.F_C.get(0, i * 3 + 0) + m_iv.F_ext.get(0, i * 3 + 0));
        m_state.a_y[i] =
            invMass * (m_iv.F_C.get(0, i * 3 + 1) + m_iv.F_ext.get(0, i * 3 + 1));
        m_state.a_theta[i] =
            invInertia * (m_iv.F_C.get(0, i * 3 + 2) + m_iv.F_ext.get(0, i * 3 + 2));
    }
}
