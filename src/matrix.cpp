#include "../include/matrix.h"

#include <assert.h>

atg_scs::Matrix::Matrix() {
    m_matrix = nullptr;
    m_data = nullptr;
    m_width = m_height = 0;
}

atg_scs::Matrix::Matrix(int width, int height, double value) {
    m_matrix = nullptr;
    m_data = nullptr;
    m_width = m_height = 0;

    initialize(width, height, value);
}

atg_scs::Matrix::~Matrix() {
    assert(m_matrix == nullptr);
}

void atg_scs::Matrix::initialize(int width, int height, double value) {
    resize(width, height);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            m_matrix[i][j] = value;
        }
    }
}

void atg_scs::Matrix::resize(int width, int height) {
    if (width == m_width && height == m_height) return;
    else if (width > m_capacityWidth || height > m_capacityHeight) {
        destroy();

        m_capacityWidth = (width > m_capacityWidth)
            ? width
            : m_capacityWidth;

        m_capacityHeight = (height > m_capacityHeight)
            ? height
            : m_capacityHeight;

        m_data = new double[(size_t)m_capacityWidth * m_capacityHeight];
        m_matrix = new double *[m_capacityHeight];
    }

    m_height = height;
    m_width = width;

    for (int i = 0; i < height; ++i) {
        m_matrix[i] = &m_data[i * width];
    }
}

void atg_scs::Matrix::destroy() {
    if (m_matrix == nullptr) {
        return;
    }

    delete[] m_matrix;
    delete[] m_data;

    m_matrix = nullptr;
    m_data = nullptr;

    m_width = m_height = 0;
}

void atg_scs::Matrix::set(int column, int row, double value) {
    m_matrix[row][column] = value;
}

void atg_scs::Matrix::set(const double *data) {
    for (int i = 0; i < m_width * m_height; ++i) {
        m_data[i] = data[i];
    }
}

double atg_scs::Matrix::get(int column, int row) {
    return m_matrix[row][column];
}

void atg_scs::Matrix::set(Matrix *reference) {
    resize(reference->m_width, reference->m_height);

    for (int i = 0; i < reference->m_height; ++i) {
        for (int j = 0; j < reference->m_width; ++j) {
            m_matrix[i][j] = reference->m_matrix[i][j];
        }
    }
}

void atg_scs::Matrix::multiply(Matrix &b, Matrix *target) {
    assert(m_width == b.m_height);

    target->resize(b.m_width, m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < b.m_width; ++j) {
            double v = 0.0;
            for (int ii = 0; ii < m_width; ++ii) {
                v += m_matrix[i][ii] * b.m_matrix[ii][j];
            }

            target->m_matrix[i][j] = v;
        }
    }
}

void atg_scs::Matrix::subtract(Matrix &b, Matrix *target) {
    assert(b.m_width == m_width);
    assert(b.m_height = m_height);

    target->resize(m_width, m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < m_width; ++j) {
            target->m_matrix[i][j] = m_matrix[i][j] - b.m_matrix[i][j];
        }
    }
}

void atg_scs::Matrix::negate(Matrix *target) {
    target->resize(m_width, m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < m_width; ++j) {
            target->m_matrix[i][j] = -m_matrix[i][j];
        }
    }
}

void atg_scs::Matrix::transpose(Matrix *target) {
    target->resize(m_height, m_width);

    for (int i = 0; i < m_width; ++i) {
        for (int j = 0; j < m_height; ++j) {
            target->m_matrix[i][j] = m_matrix[j][i];
        }
    }
}

void atg_scs::Matrix::fastRowSwap(int a, int b) {
    double *temp = m_matrix[a];
    m_matrix[a] = m_matrix[b];
    m_matrix[b] = temp;
}
