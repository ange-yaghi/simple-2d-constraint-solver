#include "../include/sle_solver.h"

atg_scs::SleSolver::SleSolver() {
    /* void */
}

atg_scs::SleSolver::~SleSolver() {
    /* void */
}

bool atg_scs::SleSolver::solve(
        Matrix &left,
        Matrix &right,
        Matrix *result,
        Matrix *previous)
{
    return false;
}

bool atg_scs::SleSolver::solveOptimized(
    Matrix &J,
    Matrix &W,
    Matrix &right,
    Matrix *result,
    Matrix *previous)
{
    return false;
}
