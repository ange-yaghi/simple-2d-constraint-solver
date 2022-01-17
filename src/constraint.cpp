#include "../include/constraint.h"

atg_scs::Constraint::Constraint(int constraintCount, int bodyCount) {
    m_constraintCount = constraintCount;
    m_bodyCount = bodyCount;

    m_index = -1; 
}

atg_scs::Constraint::~Constraint() {
    /* void */
}

void atg_scs::Constraint::calculate(Output *output, SystemState *state) {
    /* void */
}
