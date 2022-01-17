#include "../include/gaussian_elimination_sle_solver.h"

#include <cmath>
#include <assert.h>

atg_scs::GaussianEliminationSleSolver::GaussianEliminationSleSolver() {
    /* void */
}

atg_scs::GaussianEliminationSleSolver::~GaussianEliminationSleSolver() {
    /* void */
}

bool atg_scs::GaussianEliminationSleSolver::solve(
        Matrix &left,
        Matrix &right,
        Matrix *previous,
        Matrix *result)
{
    const int n = left.getWidth() + 1;
    const int m = left.getHeight();

    if (n == 0 || m == 0) return true;

    result->resize(1, m);

    Matrix A(n, m);

    for (int i = 0; i < m; ++i) {
        A.set(n - 1, i, right.get(0, i));
        for (int j = 0; j < n - 1; ++j) {
            A.set(j, i, left.get(j, i));
        }
    }

    int h = 0, k = 0;
    while (h < m && k < n) {
        int i_max = h;
        for (int i = h + 1; i < m; ++i) {
            if (A.get(k, i) > A.get(k, i_max)) {
                i_max = i;
            }
        }

        if (A.get(k, i_max) == 0) {
            ++k;
        }
        else {
            A.fastRowSwap(h, i_max);

            for (int i = h + 1; i < m; ++i) {
                const double f = A.get(k, i) / A.get(k, h);
                A.set(k, i, 0.0);

                for (int j = k + 1; j < n; ++j) {
                    A.set(j, i, A.get(j, i) - A.get(j, h) * f);
                }
            }

            ++h;
            ++k;
        }
    }

    if (A.get(n - 2, m - 1) == 0) {
        assert(false);
    }

    const double x_m = (A.get(n - 2, m - 1) == 0) 
        ? 0.0
        : A.get(n - 1, m - 1) / A.get(n - 2, m - 1);
    result->set(0, m - 1, x_m);
    for (int i = m - 2; i >= 0; --i) {
        const double b_i = A.get(n - 1, i);
        double sum = 0.0;
        for (int j = m - 1; j > i; --j) {
            sum += A.get(j, i) * result->get(0, j);
        }

        if (A.get(i, i) != 0) {
            result->set(0, i, (b_i - sum) / A.get(i, i));
        }
        else {
            result->set(0, i, 0);
        }
    }

    for (int i = 0; i < m; ++i) {
        if (std::isnan(result->get(0, i)) || std::isinf(result->get(0, i))) {
            assert(false);
        }
    }

    A.destroy();

    return true;
}
