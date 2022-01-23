#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_RIGID_BODY_SYSTEM_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_RIGID_BODY_SYSTEM_H

#include "rigid_body.h"
#include "constraint.h"
#include "force_generator.h"
#include "matrix.h"
#include "sle_solver.h"
#include "ode_solver.h"
#include "system_state.h"

#include <vector>

namespace atg_scs {
    class RigidBodySystem {
        public:
            static const int ProfilingSamples = 60 * 10;

        public:
            RigidBodySystem();
            ~RigidBodySystem();

            void initialize(SleSolver *sleSolver, OdeSolver *odeSolver);
            void reset();

            void addRigidBody(RigidBody *body);
            void removeRigidBody(RigidBody *body);
            RigidBody *getRigidBody(int i);

            void addConstraint(Constraint *constraint);
            void removeConstraint(Constraint *constraint);

            void addForceGenerator(ForceGenerator *generator);
            void removeForceGenerator(ForceGenerator *generator);

            void process(double dt, int steps = 1);

            int getRigidBodyCount() const { return (int)m_rigidBodies.size(); }
            int getConstraintCount() const { return (int)m_constraints.size(); }
            int getForceGeneratorCount() const { return (int)m_forceGenerators.size(); }

            int getFullConstraintCount() const;

            float getOdeSolveMicroseconds() const;
            float getConstraintSolveMicroseconds() const;
            float getForceEvalMicroseconds() const;
            float getConstraintEvalMicroseconds() const;

        protected:
            static float findAverage(int *samplse);

            void populateSystemState();
            void populateMassMatrices();
            void processForces();
            void processConstraints(int *evalTime, int *solveTime);

        protected:
            std::vector<RigidBody *> m_rigidBodies;
            std::vector<Constraint *> m_constraints;
            std::vector<ForceGenerator *> m_forceGenerators;

            SleSolver *m_sleSolver;
            OdeSolver *m_odeSolver;

            SystemState m_state;

            int *m_odeSolveMicroseconds;
            int *m_constraintSolveMicroseconds;
            int *m_forceEvalMicroseconds;
            int *m_constraintEvalMicroseconds;
            int m_frameIndex;

        protected:
            struct IntermediateValues {
                Matrix J, J_dot, J_T;
                Matrix M, M_inv;
                Matrix C_ks, C_kd;
                Matrix q_dot;

                Matrix reg0, reg1, reg2;

                Matrix right, left;
                Matrix F_ext, F_C;

                // Results
                Matrix lambda;
            };

            IntermediateValues m_iv;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_RIGID_BODY_SYSTEM_H */
