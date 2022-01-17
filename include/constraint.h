#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONSTRAINT_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONSTRAINT_H

#include "system_state.h"
#include "rigid_body.h"

namespace atg_scs {
    class Constraint {
        public:
            static constexpr int MaxConstraintCount = 4;
            static constexpr int MaxBodyCount = 2;

            struct Output {
                double dC_dq[MaxConstraintCount][3 * MaxBodyCount];
                double d2C_dq2[3 * MaxBodyCount][MaxConstraintCount][3 * MaxBodyCount];
                double ks[MaxConstraintCount];
                double kd[MaxConstraintCount];
                double C[MaxConstraintCount];
            };

        public:
            Constraint(int constraintCount, int bodyCount);
            virtual ~Constraint();
            
            virtual void calculate(Output *output, SystemState *state);

            int getConstraintCount() const { return m_constraintCount; }

            int m_index;

            int m_bodyCount;
            RigidBody *m_bodies[MaxBodyCount];

            double m_f_x[MaxBodyCount][MaxConstraintCount];
            double m_f_y[MaxBodyCount][MaxConstraintCount];
            double m_t[MaxBodyCount][MaxConstraintCount];

            int m_constraintCount;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONSTRAINT_H */
