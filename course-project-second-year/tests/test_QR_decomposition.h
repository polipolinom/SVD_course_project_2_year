#include <gtest/gtest.h>

#include <random>

#include "../algorithms/QR_decomposition.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "random_matrix.h"

using namespace svd_computation;

TEST(QRDecompositionTest, QRForLongDouble) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, max_number);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}

TEST(QRDecompositionTest, QRForLongDoubleMaxSize) {
    using Matrix = Matrix<long double>;

    int operations = 100;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(max_sz, max_sz, max_number);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}

TEST(QRDecompositionTest, QRForLongDoubleLowValues) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, 1.0);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}

TEST(QRDecompositionTest, QRForComplex) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, max_number);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}

TEST(QRDecompositionTest, QRForComplexMaxSize) {
    using Matrix = Matrix<Complex>;

    int operations = 100;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(max_sz, max_sz, 1.0);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}

TEST(QRDecompositionTest, QRForComplexLowValues) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, 1.0);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}
