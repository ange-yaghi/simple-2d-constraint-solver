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
    destroy();

    m_height = height;
    m_width = width;

    m_data = new double[(size_t)width * height];
    m_matrix = new double *[height];
    for (int i = 0; i < height; ++i) {
        m_matrix[i] = &m_data[i * width];
    }

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            m_matrix[i][j] = value;
        }
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
    assert(m_width == reference->m_width);
    assert(m_height == reference->m_height);

    for (int i = 0; i < reference->m_height; ++i) {
        for (int j = 0; j < reference->m_width; ++j) {
            m_matrix[i][j] = reference->m_matrix[i][j];
        }
    }
}

void atg_scs::Matrix::multiply(Matrix &b, Matrix *target) {
    assert(m_width == b.m_height);
    assert(target->m_width == b.m_width);
    assert(target->m_height == m_height);

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
    assert(target->m_width == m_width);
    assert(target->m_height == m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < m_width; ++j) {
            target->m_matrix[i][j] = m_matrix[i][j] - b.m_matrix[i][j];
        }
    }
}

void atg_scs::Matrix::negate(Matrix *target) {
    assert(target->m_width == m_width);
    assert(target->m_height == m_height);

    for (int i = 0; i < m_height; ++i) {
        for (int j = 0; j < m_width; ++j) {
            target->m_matrix[i][j] = -m_matrix[i][j];
        }
    }
}

void atg_scs::Matrix::transpose(Matrix *target) {
    assert(target->m_width == m_height);
    assert(target->m_height == m_width);

    for (int i = 0; i < m_width; ++i) {
        for (int j = 0; j < m_height; ++j) {
            target->m_matrix[i][j] = m_matrix[j][i];
        }
    }
}
