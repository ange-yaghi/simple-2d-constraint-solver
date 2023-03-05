#ifndef ATG_SIMPLE_2D_CONSTRAINT_SIMPLE_GEAR_CONSTRAINT_H
#define ATG_SIMPLE_2D_CONSTRAINT_SIMPLE_GEAR_CONSTRAINT_H

#include "constraint.h"

namespace atg_scs {
    class SimpleGearConstraint : public Constraint {
        public:
            SimpleGearConstraint();
            virtual ~SimpleGearConstraint();

            void setBody1(RigidBody *body) { m_bodies[0] = body; }
            void setBody2(RigidBody *body) { m_bodies[1] = body; }

            virtual void calculate(Output *output, SystemState *system);

            double m_ratio;
            bool m_neutral;

            double m_ks;
            double m_kd;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SIMPLE_GEAR_CONSTRAINT_H */
