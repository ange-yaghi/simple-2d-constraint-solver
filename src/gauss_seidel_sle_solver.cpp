#include "../include/gauss_seidel_sle_solver.h"

#include <cmath>
#include <assert.h>

atg_scs::GaussSeidelSleSolver::GaussSeidelSleSolver() {
    m_maxIterations = 1000;
    m_minDelta = 1E-15;
}

atg_scs::GaussSeidelSleSolver::~GaussSeidelSleSolver() {
    /* void */
}

bool atg_scs::GaussSeidelSleSolver::solve(
        Matrix &left,
        Matrix &right,
        Matrix *previous,
        Matrix *result)
{
    const int n = right.getHeight();
    
    result->resize(1, n);

    if (previous != nullptr && previous->getHeight() == n) {
        result->set(previous);
    }

    for (int i = 0; i < m_maxIterations; ++i) {
        const double maxDelta = solveIteration(
                left,
                right,
                result,
                result);

        if (maxDelta < m_minDelta) {
            return true;
        }
    }

    return false;
}

double atg_scs::GaussSeidelSleSolver::solveIteration(
        Matrix &left,
        Matrix &right,
        Matrix *k_next,
        Matrix *k)
{
    double maxDifference = 0.0;
    const int n = k->getHeight();

    for (int i = 0; i < n; ++i) {
        double s0 = 0.0, s1 = 0.0;
        for (int j = 0; j < i; ++j) {
            s0 += left.get(j, i) * k_next->get(0, j);
        }

        for (int j = i + 1; j < n; ++j) {
            s1 += left.get(j, i) * k->get(0, j);
        }

        const double k_next_i =
            (1 / left.get(i, i)) * (right.get(0, i) - s0 - s1);
        const double delta = std::abs(k->get(0, i) - k_next_i);
        maxDifference = std::fmax(delta, maxDifference);

        k_next->set(0, i, k_next_i);
    }

    return maxDifference;
}
