#include "../include/sparse_matrix.h"

#include "../include/matrix.h"

#include <algorithm>
#include <assert.h>

atg_scs::SparseMatrix::SparseMatrix() {
    m_matrix = nullptr;
    m_data = nullptr;
    m_blockData = nullptr;
    m_width = m_height = 0;
    m_capacityHeight = 0;
    m_stride = 1;
    m_elementsPerRow = 0;
}

atg_scs::SparseMatrix::~SparseMatrix() {
    assert(m_matrix == nullptr);
    assert(m_data == nullptr);
    assert(m_blockData == nullptr);
}

void atg_scs::SparseMatrix::initialize(
        int width,
        int height,
        int stride,
        int entries)
{
    m_stride = stride;
    m_elementsPerRow = entries;

    resize(width, height);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < m_elementsPerRow; ++j) {
            setEmpty(i, j);
        }
    }
}

void atg_scs::SparseMatrix::resize(int width, int height) {
    if (width == m_width && height == m_height) return;
    else if (height > m_capacityHeight) {
        destroy();

        m_capacityHeight = (height > m_capacityHeight)
            ? height
            : m_capacityHeight;

        m_data = new double[(size_t)m_stride * m_elementsPerRow * m_capacityHeight];
        m_matrix = new double *[m_capacityHeight];
        m_blockData = new uint8_t[(size_t)m_capacityHeight * m_elementsPerRow];
    }

    m_height = height;
    m_width = width;

    for (int i = 0; i < height; ++i) {
        m_matrix[i] = &m_data[i * m_elementsPerRow * m_stride];
    }
}

void atg_scs::SparseMatrix::destroy() {
    if (m_matrix == nullptr) {
        return;
    }

    delete[] m_matrix;
    delete[] m_data;
    delete[] m_blockData;

    m_matrix = nullptr;
    m_data = nullptr;
    m_blockData = nullptr;

    m_width = m_height = 0;
}

void atg_scs::SparseMatrix::expand(Matrix *matrix) {
    matrix->initialize(m_width, m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < m_elementsPerRow; ++j) {
            const uint8_t block = m_blockData[i * m_elementsPerRow + j];
            if (block == 0xFF) continue;
            else {
                for (int k = 0; k < m_stride; ++k) {
                    matrix->set(block, i, m_matrix[i][j * m_stride + k]);
                }
            }
        }
    }
}

void atg_scs::SparseMatrix::multiplyTranspose(SparseMatrix &b_T, Matrix *target) {
    assert(m_width == b_T.m_width);
    assert(m_stride == b_T.m_stride);

    target->initialize(b_T.m_height, m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < b_T.m_height; ++j) {
            double dot = 0;
            for (int k = 0; k < m_elementsPerRow; ++k) {
                const uint8_t block0 = m_blockData[i * m_elementsPerRow + k];
                for (int l = 0; l < m_elementsPerRow; ++l) {
                    const uint8_t block1 = b_T.m_blockData[j * b_T.m_elementsPerRow + l];
                    if (block0 == block1) {
                        for (int m = 0; m < m_stride; ++m) {
                            dot +=
                                m_matrix[i][k * m_stride + m]
                                * b_T.m_matrix[j][l * m_stride + m];
                        }
                    }
                }
            }

            target->set(j, i, dot);
        }
    }
}

void atg_scs::SparseMatrix::rightScale(Matrix &scale, SparseMatrix *target) {
    assert(scale.getWidth() == 1);
    assert(scale.getHeight() == m_width);

    target->initialize(m_width, m_height, m_stride, m_elementsPerRow);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < m_elementsPerRow; ++j) {
            target->m_matrix[i][j] = scale.get(0, j) * m_matrix[i][j];
        }
    }
}
