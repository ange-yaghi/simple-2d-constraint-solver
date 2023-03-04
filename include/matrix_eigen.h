#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_MATRIX_EIGEN_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_MATRIX_EIGEN_H

#include <assert.h>

#ifdef ATG_S2C_USE_EIGEN

#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT

#include "Eigen/Dense"

namespace atg_scs
{
    class Matrix
    {
        typedef Eigen::MatrixXd MatrixType;

    public:
        Matrix();
        Matrix(int width, int height, double value = 0.0);
        ~Matrix();

        void initialize(int width, int height, double value);
        void initialize(int width, int height);
        void resize(int width, int height);
        void destroy();

        void set(const double *data);

        __forceinline void set(int column, int row, double value)
        {
            m_matrix(row, column) = value;
        }

        __forceinline void add(int column, int row, double value)
        {
            m_matrix(row, column) += value;
        }

        __forceinline double get(int column, int row)
        {
            return m_matrix(row, column);
        }

        void set(Matrix *reference);

        void multiply(Matrix &b, Matrix *target);
        void componentMultiply(Matrix &b, Matrix *target);
        void transposeMultiply(Matrix &b, Matrix *target);
        void leftScale(Matrix &scale, Matrix *target);
        void rightScale(Matrix &scale, Matrix *target);
        void scale(double s, Matrix *target);
        void subtract(Matrix &b, Matrix *target);
        void add(Matrix &b, Matrix *target);
        void negate(Matrix *target);
        bool equals(Matrix &b, double err = 1e-6);
        double vectorMagnitudeSquared() const;
        double dot(Matrix &b) const;

        void madd(Matrix &b, double s);
        void pmadd(Matrix &b, double s);

        void transpose(Matrix *target);
        int getWidth() const { return m_matrix.cols(); }
        int getHeight() const { return m_matrix.rows(); }

        __forceinline void fastRowSwap(int a, int b)
        {
            auto row_a = m_matrix.row(a);
            auto row_b = m_matrix.row(a);
            m_matrix.row(a) = row_b;
            m_matrix.row(b) = row_a;
        }

    protected:
        MatrixType m_matrix;
    };
} /* namespace atg_scs */

#endif

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_MATRIX_H */
