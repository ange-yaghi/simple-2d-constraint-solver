#include "../include/matrix_eigen.h"

#ifdef ATG_S2C_USE_EIGEN
#include <algorithm>
#include <cassert>
#include <iostream>

atg_scs::Matrix::Matrix()
{
}

atg_scs::Matrix::Matrix(int width, int height, double value)
{
    m_matrix = MatrixType::Constant(height, width, value);
}

atg_scs::Matrix::~Matrix()
{
}

void atg_scs::Matrix::initialize(int width, int height, double value)
{
    m_matrix = MatrixType::Constant(height, width, value);
}

void atg_scs::Matrix::initialize(int width, int height)
{
    m_matrix = MatrixType::Zero(height, width);
}

void atg_scs::Matrix::resize(int width, int height)
{
    m_matrix = MatrixType::Zero(height, width);
}

void atg_scs::Matrix::destroy()
{
}

void atg_scs::Matrix::set(const double *data)
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MaxtrixRowMajor;

    m_matrix = Eigen::Map<const MaxtrixRowMajor>(data, m_matrix.rows(), m_matrix.cols());
}

void atg_scs::Matrix::set(Matrix *reference)
{
    m_matrix = reference->m_matrix;
}

void atg_scs::Matrix::multiply(Matrix &b, Matrix *target)
{
    target->m_matrix = m_matrix * b.m_matrix;
}

void atg_scs::Matrix::componentMultiply(Matrix &b, Matrix *target)
{
    target->m_matrix = (b.m_matrix.array() * m_matrix.array()).matrix();
}

void atg_scs::Matrix::transposeMultiply(Matrix &b, Matrix *target)
{
    target->m_matrix = m_matrix.transpose() * b.m_matrix;
}

void atg_scs::Matrix::leftScale(Matrix &scale, Matrix *target)
{
    // Scale is a column vector (scale.cols == 1 && scale.rows() == rows())
    assert(scale.getWidth() == 1);
    assert(scale.getHeight() == getHeight());

    auto vector = scale.m_matrix.col(0).array();
    target->m_matrix = (m_matrix.array().colwise() * vector).matrix();
}

void atg_scs::Matrix::rightScale(Matrix &scale, Matrix *target)
{
    assert(scale.getWidth() == 1);
    assert(scale.getHeight() == getWidth());

    auto vector = scale.m_matrix.col(0).transpose().array();
    target->m_matrix = (m_matrix.array().rowwise() * vector).matrix();
}

void atg_scs::Matrix::scale(double s, Matrix *target)
{
    target->m_matrix = m_matrix * s;
}

void atg_scs::Matrix::subtract(Matrix &b, Matrix *target)
{
    target->m_matrix = m_matrix - b.m_matrix;
}

void atg_scs::Matrix::add(Matrix &b, Matrix *target)
{
    target->m_matrix = m_matrix + b.m_matrix;
}

void atg_scs::Matrix::negate(Matrix *target)
{
    target->m_matrix = -m_matrix;
}

bool atg_scs::Matrix::equals(Matrix &b, double err)
{
    if (getWidth() != b.getWidth())
        return false;

    if (getHeight() != b.getHeight())
        return false;

    return ((b.m_matrix - m_matrix).array() < err).all();
}

double atg_scs::Matrix::vectorMagnitudeSquared() const
{
    assert(getWidth() == 1);
    return m_matrix.col(0).squaredNorm();
}

double atg_scs::Matrix::dot(Matrix &b) const
{
    assert(getWidth() == 1);
    assert(b.getWidth() == 1);
    assert(b.getHeight() == getHeight());

    return m_matrix.col(0).dot(b.m_matrix.col(0));
    ;
}

void atg_scs::Matrix::madd(Matrix &b, double s)
{
    m_matrix += b.m_matrix * s;
}

void atg_scs::Matrix::pmadd(Matrix &b, double s)
{
    m_matrix = s * m_matrix + b.m_matrix;
}

void atg_scs::Matrix::transpose(Matrix *target)
{
    target->m_matrix = m_matrix.transpose();
}

#endif