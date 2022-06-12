#include "../include/conjugate_gradient_sle_solver.h"

#include <cmath>
#include <assert.h>

atg_scs::ConjugateGradientSleSolver::ConjugateGradientSleSolver()
    : atg_scs::SleSolver(false)
{
    m_maxIterations = 1000;
    m_maxError = 1E-2;
    m_minError = 1E-3;
}

atg_scs::ConjugateGradientSleSolver::~ConjugateGradientSleSolver() {
    m_mreg0.destroy();
    m_mreg1.destroy();
    m_Ap.destroy();
    m_x.destroy();
    m_r.destroy();
    m_p.destroy();
}

bool atg_scs::ConjugateGradientSleSolver::solve(
        SparseMatrix<3> &J,
        Matrix &W,
        Matrix &right,
        Matrix *previous,
        Matrix *result)
{
    const int n = right.getHeight();

    m_r.resize(1, n);
    m_p.resize(1, n);
    m_Ap.resize(1, n);
    m_x.initialize(1, n);

    if (previous != nullptr && previous->getHeight() == n) {
        m_x.set(previous);
    }

    result->resize(1, n);

    m_r.set(&right);
    multiply(J, W, m_x, &m_Ap);
    m_r.madd(m_Ap, -1);

    if (sufficientlySmall(m_r, right)) {
        goto succeeded;
    }

    m_p.set(&m_r);

    for (int k = 0; k < m_maxIterations; ++k) {
        multiply(J, W, m_p, &m_Ap);

        const double rk_mag = m_r.vectorMagnitudeSquared();
        const double alpha = rk_mag / m_p.dot(m_Ap);
        m_x.madd(m_p, alpha);
        m_r.madd(m_Ap, -alpha);

        if (sufficientlySmall(m_r, right)) {
            goto succeeded;
        }

        const double rk1_mag = m_r.vectorMagnitudeSquared();
        const double beta = rk1_mag / rk_mag;
        m_p.pmadd(m_r, beta);
    }

    return false;

succeeded:
    result->set(&m_x);

    return true;
}

void atg_scs::ConjugateGradientSleSolver::multiply(
    SparseMatrix<3> &J,
    Matrix &W,
    Matrix &x,
    Matrix *target)
{
    // A = J * W * J_T
    target->resize(1, x.getHeight());
    J.transposeMultiplyVector(x, &m_mreg0);
    W.componentMultiply(m_mreg0, &m_mreg1);
    J.multiply(m_mreg1, target);
}

bool atg_scs::ConjugateGradientSleSolver::sufficientlySmall(
    Matrix &x,
    Matrix &target) const
{
    for (int i = 0; i < x.getHeight(); ++i) {
        const double err = x.get(0, i);
        const double t = target.get(0, i);
        if (std::abs(err) > std::fmax(std::abs(m_maxError * t), m_minError)) {
            return false;
        }
    }

    return true;
}
