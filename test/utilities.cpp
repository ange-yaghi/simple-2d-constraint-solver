#include "utilities.h"

#include <gtest/gtest.h>

void compareMatrix(atg_scs::Matrix &a, atg_scs::Matrix &b, double err) {
    EXPECT_EQ(a.getWidth(), b.getWidth());
    EXPECT_EQ(a.getHeight(), b.getHeight());

    for (int i = 0; i < a.getWidth(); ++i) {
        for (int j = 0; j < b.getWidth(); ++j) {
            ASSERT_NEAR(a.get(i, j), b.get(i, j), err);
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
