#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_GAUSSIAN_ELIMINATION_SLE_SOLVER_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_GAUSSIAN_ELIMINATION_SLE_SOLVER_H

#include "sle_solver.h"

namespace atg_scs {
    class GaussianEliminationSleSolver : public SleSolver {
        public:
            GaussianEliminationSleSolver();
            virtual ~GaussianEliminationSleSolver();

            virtual bool solve(
                    SparseMatrix<3> &J,
                    Matrix &W,
                    Matrix &right,
                    Matrix *result,
                    Matrix *previous);

        protected:
            Matrix m_a;
            Matrix m_M;
            SparseMatrix<3> m_reg;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_GAUSSIAN_ELIMINATION_SLE_SOLVER_H */
