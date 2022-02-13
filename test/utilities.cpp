#include "utilities.h"

#include <gtest/gtest.h>

void compareMatrix(atg_scs::Matrix &a, atg_scs::Matrix &b, double err) {
    EXPECT_EQ(a.getWidth(), b.getWidth());
    EXPECT_EQ(a.getHeight(), b.getHeight());

    for (int i = 0; i < a.getHeight(); ++i) {
        for (int j = 0; j < b.getWidth(); ++j) {
            ASSERT_NEAR(a.get(j, i), b.get(j, i), err);
        }
    }
}

void fullToSparse(atg_scs::Matrix &full, atg_scs::SparseMatrix<3, 2> *target, int stride) {
    const int entries = full.getWidth() / stride;
    target->initialize(full.getWidth(), full.getHeight());

    for (int i = 0; i < full.getHeight(); ++i) {
        int currentEntry = 0;
        for (int j = 0; j < entries; ++j) {
            bool nonzero = false;
            for (int k = 0; k < stride; ++k) {
                if (full.get(j * stride + k, i) != 0) {
                    nonzero = true;
                    break;
                }
            }

            if (nonzero) {
                const int entry = currentEntry++;

                for (int slice = 0; slice < stride; ++slice) {
                    target->setBlock(i, entry, j);
                    target->set(i, entry, slice, full.get(j * stride + slice, i));
                }
            }
        }
    }
}

void JWJ_t(atg_scs::Matrix &J, atg_scs::Matrix &W, atg_scs::Matrix *target) {
    atg_scs::Matrix temp;
    atg_scs::Matrix J_T;

    J.transpose(&J_T);

    J.rightScale(W, &temp);
    temp.multiply(J_T, target);

    temp.destroy();
    J_T.destroy();
}
