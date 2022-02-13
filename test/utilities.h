#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_UTILITIES_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_UTILITIES_H

#include "../include/matrix.h"
#include "../include/sparse_matrix.h"

void compareMatrix(atg_scs::Matrix &a, atg_scs::Matrix &b, double err = 1E-6);
void fullToSparse(atg_scs::Matrix &full, atg_scs::SparseMatrix<3, 2> *target, int stride);
void JWJ_t(atg_scs::Matrix &J, atg_scs::Matrix &W, atg_scs::Matrix *target);

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_UTILITIES_H */
