#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H

#include <assert.h>
#include <stdint.h>

namespace atg_scs {
    class Matrix;
    class SparseMatrix {
        public:
            SparseMatrix();
            ~SparseMatrix();

            void initialize(
                    int width,
                    int height,
                    int stride,
                    int entries);
            void resize(int width, int height);
            void destroy();

            void expand(Matrix *matrix);

            inline void setBlock(int row, int entry, uint8_t index) {
                assert(row >= 0 && row < m_height);
                assert(entry >= 0 && entry < m_elementsPerRow);
                assert(index < m_width);

                m_blockData[row * m_elementsPerRow + entry] = index;
            }

            inline void set(int row, int entry, int slice, double v) {
                assert(row >= 0 && row < m_height);
                assert(entry >= 0 && entry < m_elementsPerRow);
                assert(slice < m_stride);

                m_matrix[row][entry * m_stride + slice] = v;
            }

            inline void setEmpty(int row, int col) {
                assert(row >= 0 && row < m_height);
                assert(col >= 0 && col < m_elementsPerRow);

                m_blockData[row * m_elementsPerRow + col] = 0xFF;
                for (int i = 0; i < m_stride; ++i) {
                    m_matrix[row][col * m_stride + i] = 0;
                }
            }

            void multiplyTranspose(SparseMatrix &b_T, Matrix *target);
            void rightScale(Matrix &scale, SparseMatrix *target);

            int getWidth() const { return m_width; }
            int getHeight() const { return m_height; }

        protected:
            double **m_matrix;
            double *m_data;
            uint8_t *m_blockData;

            int m_stride;
            int m_elementsPerRow;
            int m_width;
            int m_height;
            int m_capacityHeight;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H */
