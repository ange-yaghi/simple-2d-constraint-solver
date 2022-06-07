#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CLUTCH_CONSTRAINT_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CLUTCH_CONSTRAINT_H

#include "constraint.h"

namespace atg_scs {
    class CluchConstraint : public Constraint {
        public:
            CluchConstraint();
            virtual ~CluchConstraint();
            
            void setBody1(RigidBody *body) { m_bodies[0] = body; }
            void setBody2(RigidBody *body) { m_bodies[1] = body; }

            virtual void calculate(Output *output, SystemState *system);
            virtual void limit(atg_scs::Matrix *lambda, atg_scs::SystemState *state);

            double m_maxTorque;
            double m_minTorque;

            double m_ks;
            double m_kd;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CLUTCH_CONSTRAINT_H */
