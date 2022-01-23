#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONSTRAINT_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONSTRAINT_H

#include "system_state.h"
#include "rigid_body.h"

namespace atg_scs {
    class Constraint {
        public:
            static constexpr int MaxConstraintCount = 6;
            static constexpr int MaxBodyCount = 2;

            struct Output {
                double C[MaxConstraintCount];
                double J[MaxConstraintCount][3 * MaxBodyCount];
                double J_dot[MaxConstraintCount][3 * MaxBodyCount];
                double ks[MaxConstraintCount];
                double kd[MaxConstraintCount];
            };

            struct Lambda {
                double lambda[MaxConstraintCount];
            };

        public:
            Constraint(int constraintCount, int bodyCount);
            virtual ~Constraint();

            virtual void calculate(Output *output, SystemState *state);
            virtual int getConstraintCount() const { return m_constraintCount; }

            int m_index;
            int m_bodyCount;
            RigidBody *m_bodies[MaxBodyCount];

        protected:
            int m_constraintCount;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONSTRAINT_H */
