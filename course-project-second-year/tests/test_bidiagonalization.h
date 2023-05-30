#include <gtest/gtest.h>

#include <random>

#include "../algorithms/bidiagonalization.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "random_matrix.h"

using namespace svd_computation;

TEST(BidiagonalizationTest, BidiagonalizeLongDouble) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B));
        EXPECT_TRUE(details::is_unitary(left));
        EXPECT_TRUE(details::is_unitary(right));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A));
    }
}

TEST(BidiagonalizationTest, BidiagonalizeLongDoubleMaxSize) {
    using Matrix = Matrix<long double>;

    int operations = 100;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(max_sz, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B));
        EXPECT_TRUE(details::is_unitary(left));
        EXPECT_TRUE(details::is_unitary(right));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A));
    }
}

TEST(BidiagonalizationTest, BidiagonalizeLongDoubleLowValues) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, 1.0);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B));
        EXPECT_TRUE(details::is_unitary(left));
        EXPECT_TRUE(details::is_unitary(right));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A));
    }
}

TEST(BidiagonalizationTest, bidiagonalizeComplex) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B));
        EXPECT_TRUE(details::is_unitary(left));
        EXPECT_TRUE(details::is_unitary(right));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A));

        // check B is real matrix
        for (size_t i = 0; i < B.height(); ++i) {
            for (size_t j = 0; j < B.width(); ++j) {
                EXPECT_NEAR(B(i, j).Im(), 0, constants::DEFAULT_EPSILON);
            }
        }
    }
}

TEST(BidiagonalizationTest, bidiagonalizeComplexMaxSize) {
    using Matrix = Matrix<Complex>;

    int operations = 100;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(max_sz, max_sz, max_number);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B));
        EXPECT_TRUE(details::is_unitary(left));
        EXPECT_TRUE(details::is_unitary(right));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A));

        // check B is real matrix
        for (size_t i = 0; i < B.height(); ++i) {
            for (size_t j = 0; j < B.width(); ++j) {
                EXPECT_NEAR(B(i, j).Im(), 0, constants::DEFAULT_EPSILON);
            }
        }
    }
}

TEST(BidiagonalizationTest, bidiagonalizeLowValues) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, 1.0);

        Matrix left, right;
        auto B = bidiagonalize(A, &left, &right);

        EXPECT_TRUE(details::is_bidiagonal(B));
        EXPECT_TRUE(details::is_unitary(left));
        EXPECT_TRUE(details::is_unitary(right));
        EXPECT_TRUE(details::is_zero(left * B * conjugate(right) - A));

        // check B is real matrix
        for (size_t i = 0; i < B.height(); ++i) {
            for (size_t j = 0; j < B.width(); ++j) {
                EXPECT_NEAR(B(i, j).Im(), 0, constants::DEFAULT_EPSILON);
            }
        }
    }
}
