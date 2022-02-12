#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_UTILITIES_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_UTILITIES_H

#include "../include/matrix.h"

void compareMatrix(atg_scs::Matrix &a, atg_scs::Matrix &b, double err = 1E-6);
void JWJ_t(atg_scs::Matrix &J, atg_scs::Matrix &W, atg_scs::Matrix *target);

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_UTILITIES_H */
