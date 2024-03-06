#pragma once

#include <gtest/gtest.h>

#include "../algorithms/QR_decomposition.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "random_objects.h"

using namespace svd_computation;

TEST(QRDecompositionTest, QRForLongDouble) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, max_number);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R, eps));
        EXPECT_TRUE(details::is_unitary(Q, eps));
        EXPECT_TRUE(details::is_zero(Q * R - A, eps));
    }
}

TEST(QRDecompositionTest, QRForLongDoubleMaxSize) {
    using Matrix = Matrix<long double>;

    int operations = 500;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(max_sz, max_sz, max_number);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R, eps));
        EXPECT_TRUE(details::is_unitary(Q, eps));
        EXPECT_TRUE(details::is_zero(Q * R - A, eps));
    }
}

TEST(QRDecompositionTest, QRForLongDoubleLowValues) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, 1.0);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R, eps));
        EXPECT_TRUE(details::is_unitary(Q, eps));
        EXPECT_TRUE(details::is_zero(Q * R - A, eps));
    }
}

TEST(QRDecompositionTest, QRForComplex) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, max_number);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R, eps));
        EXPECT_TRUE(details::is_unitary(Q, eps));
        EXPECT_TRUE(details::is_zero(Q * R - A, eps));
    }
}

TEST(QRDecompositionTest, QRForComplexMaxSize) {
    using Matrix = Matrix<Complex>;

    int operations = 500;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(max_sz, max_sz, 1.0);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R, eps));
        EXPECT_TRUE(details::is_unitary(Q, eps));
        EXPECT_TRUE(details::is_zero(Q * R - A, eps));
    }
}

TEST(QRDecompositionTest, QRForComplexLowValues) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, 1.0);

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R, eps));
        EXPECT_TRUE(details::is_unitary(Q, eps));
        EXPECT_TRUE(details::is_zero(Q * R - A, eps));
    }
}
