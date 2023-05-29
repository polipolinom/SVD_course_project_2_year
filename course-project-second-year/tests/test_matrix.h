#include <gtest/gtest.h>

#include "../types/complex.h"
#include "../types/matrix.h"
#include "../types/vector.h"

using namespace svd_computation;

TEST(MatrixTest, MatrixMult) {
    Matrix<long double> A = {{-4.5, -3.6, 7.8}, {5.2, -8.9, 1.6}, {-0.3, 4, 5}};
    Matrix<long double> B = {{0.9, 4.67, -8.95}, {3, -4.5, 234}, {1.23, 5.4, 3}};

    Matrix<long double> AB = {{-5.256, 37.305, -778.725}, {-20.052, 72.974, -2124.34}, {17.88, 7.599, 953.685}};
    Matrix<long double> BA = {{22.919, -80.603, -30.258}, {-107.1, 965.25, 1186.2}, {21.645, -40.488, 33.234}};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ((A * B)(i, j), AB(i, j));
            EXPECT_DOUBLE_EQ((B * A)(i, j), BA(i, j));
        }
    }
}

TEST(MatrixTest, VectorMatrixMult) {
    Matrix<long double> A = {{-4.5, -3.6, 7.8}, {5.2, -8.9, 1.6}, {-0.3, 4, 5}};
    Vector<long double> v = {1.2, 3.4, 1.234};
    Vector<long double> u({234, 6.78, 0.98}, Vector<long double>::Orientation::Horizontal);

    Vector<long double> Av = {-8.0148, -22.0456, 19.41};
    Vector<long double> uA({-1018.038, -898.822, 1840.948}, Vector<long double>::Orientation::Horizontal);

    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ((A * v)[i], Av[i]);
        EXPECT_DOUBLE_EQ((u * A)[i], uA[i]);
    }
}

TEST(MatrixTest, Transpose) {
    Matrix<int> A = {{1, 2, 3}, {4, 5, 6}};
    EXPECT_EQ(A.height(), 2);
    EXPECT_EQ(A.width(), 3);

    auto B = transpose(A);

    EXPECT_EQ(B.height(), 3);
    EXPECT_EQ(B.width(), 2);

    EXPECT_EQ(A.height(), 2);
    EXPECT_EQ(A.width(), 3);

    A.transpose();

    EXPECT_EQ(A.height(), 3);
    EXPECT_EQ(A.width(), 2);

    A.transpose();

    EXPECT_EQ(A.height(), 2);
    EXPECT_EQ(A.width(), 3);
}