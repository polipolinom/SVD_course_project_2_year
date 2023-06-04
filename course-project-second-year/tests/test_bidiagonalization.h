#pragma once

#include <gtest/gtest.h>

#include "../algorithms/bidiagonalization.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "random_objects.h"

using namespace svd_computation;

TEST(BidiagonalizationTest, BidiagonalizeLongDouble) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 300;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B, eps));
        EXPECT_TRUE(details::is_unitary(left, eps));
        EXPECT_TRUE(details::is_unitary(right, eps));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A, eps));
    }
}

TEST(BidiagonalizationTest, BidiagonalizeLongDoubleMaxSize) {
    using Matrix = Matrix<long double>;

    int operations = 100;
    int max_sz = 300;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(max_sz, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B, eps));
        EXPECT_TRUE(details::is_unitary(left, eps));
        EXPECT_TRUE(details::is_unitary(right, eps));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A, eps));
    }
}

TEST(BidiagonalizationTest, BidiagonalizeLongDoubleLowValues) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 300;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, 1.0);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B, eps));
        EXPECT_TRUE(details::is_unitary(left, eps));
        EXPECT_TRUE(details::is_unitary(right, eps));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A, eps));
    }
}

TEST(BidiagonalizationTest, bidiagonalizeComplex) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 300;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B, eps));
        EXPECT_TRUE(details::is_unitary(left, eps));
        EXPECT_TRUE(details::is_unitary(right, eps));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A, eps));

        // check B is real matrix
        for (size_t i = 0; i < B.height(); ++i) {
            for (size_t j = 0; j < B.width(); ++j) {
                EXPECT_NEAR(B(i, j).Im(), 0, eps);
            }
        }
    }
}

TEST(BidiagonalizationTest, bidiagonalizeComplexMaxSize) {
    using Matrix = Matrix<Complex>;

    int operations = 100;
    int max_sz = 300;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(max_sz, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B, eps));
        EXPECT_TRUE(details::is_unitary(left, eps));
        EXPECT_TRUE(details::is_unitary(right, eps));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A, eps));

        // check B is real matrix
        for (size_t i = 0; i < B.height(); ++i) {
            for (size_t j = 0; j < B.width(); ++j) {
                EXPECT_NEAR(B(i, j).Im(), 0, eps);
            }
        }
    }
}

TEST(BidiagonalizationTest, bidiagonalizeComplexLowValues) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 300;
    long double max_number = 1e5;

    long double eps = 1e-10;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, 1.0);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B, eps));
        EXPECT_TRUE(details::is_unitary(left, eps));
        EXPECT_TRUE(details::is_unitary(right, eps));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A, eps));

        // check B is real matrix
        for (size_t i = 0; i < B.height(); ++i) {
            for (size_t j = 0; j < B.width(); ++j) {
                EXPECT_NEAR(B(i, j).Im(), 0, eps);
            }
        }
    }
}
