#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SEIDEL_SLE_SOLVER_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SEIDEL_SLE_SOLVER_H

#include "sle_solver.h"

namespace atg_scs {
    class GaussSeidelSleSolver : public SleSolver {
        public:
            GaussSeidelSleSolver();
            virtual ~GaussSeidelSleSolver();

            virtual bool solve(
                    Matrix &left,
                    Matrix &right,
                    Matrix *result,
                    Matrix *previous);

            int m_maxIterations;
            double m_minDelta;

        protected:
            double solveIteration(
                    Matrix &left,
                    Matrix &right,
                    Matrix *result,
                    Matrix *previous);
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SEIDEL_SLE_SOLVER_H */
