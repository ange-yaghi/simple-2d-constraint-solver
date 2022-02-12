#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SLE_SOLVER_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SLE_SOLVER_H

#include "matrix.h"

namespace atg_scs {
    class SleSolver {
        public:
            SleSolver();
            virtual ~SleSolver();

            virtual bool solve(
                    Matrix &J,
                    Matrix &W,
                    Matrix &right,
                    Matrix *result,
                    Matrix *previous);
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SLE_SOLVER_H */
