#include "../include/gaussian_elimination_sle_solver.h"

#include <cmath>
#include <assert.h>
#include <fstream>

atg_scs::GaussianEliminationSleSolver::GaussianEliminationSleSolver() {
    m_a.initialize(1, 1);
}

atg_scs::GaussianEliminationSleSolver::~GaussianEliminationSleSolver() {
    m_a.destroy();
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

    Matrix &A = m_a;
    A.resize(n, m);

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
            if (std::abs(A.get(k, i)) > std::abs(A.get(k, i_max))) {
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
        std::fstream f("bad_matrix.csv", std::ios::out);
        for (int i = 0; i < left.getHeight(); ++i) {
            for (int j = 0; j < left.getWidth(); ++j) {
                f << left.get(j, i);
                if (j < left.getWidth() - 1) f << ",";
                else f << "\n";
            }
        }

        f.close();
        f.open("rref.csv", std::ios::out);
        for (int i = 0; i < A.getHeight(); ++i) {
            for (int j = 0; j < A.getWidth(); ++j) {
                f << A.get(j, i);
                if (j < A.getWidth() - 1) f << ",";
                else f << "\n";
            }
        }

        f.close();
        f.open("bad_vector.csv", std::ios::out);
        for (int i = 0; i < right.getHeight(); ++i) {
            for (int j = 0; j < right.getWidth(); ++j) {
                f << right.get(j, i);
                if (j < right.getWidth() - 1) f << ",";
                else f << "\n";
            }
        }

        assert(false);
    }

    const double x_m = A.get(n - 1, m - 1) / A.get(n - 2, m - 1);
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
            std::fstream f("bad_matrix.csv", std::ios::out);
            for (int i = 0; i < left.getHeight(); ++i) {
                for (int j = 0; j < left.getWidth(); ++j) {
                    f << left.get(j, i);
                    if (j < left.getWidth() - 1) f << ",";
                    else f << "\n";
                }
            }

            f.close();
            f.open("rref.csv", std::ios::out);
            for (int i = 0; i < A.getHeight(); ++i) {
                for (int j = 0; j < A.getWidth(); ++j) {
                    f << A.get(j, i);
                    if (j < A.getWidth() - 1) f << ",";
                    else f << "\n";
                }
            }

            f.close();
            f.open("bad_vector.csv", std::ios::out);
            for (int i = 0; i < right.getHeight(); ++i) {
                for (int j = 0; j < right.getWidth(); ++j) {
                    f << right.get(j, i);
                    if (j < right.getWidth() - 1) f << ",";
                    else f << "\n";
                }
            }
            f.close();

            assert(false);
        }
    }

    return true;
}
