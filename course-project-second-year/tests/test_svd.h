#pragma once

#include <gtest/gtest.h>

#include "../algorithms/svd_computation.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "Eigen/eigen"
#include "random_objects.h"

using namespace svd_computation;

TEST(SVDTest, TestForLongDouble) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1, max_sz, max_number);

        Matrix left, right;
        auto sigma = compute_svd(A, &left, &right);

        // check in sigma values > 0
        for (auto s : sigma) {
            EXPECT_TRUE(s >= 0);
        }

        int mx = std::max(A.height(), A.width());

        EXPECT_TRUE(details::is_unitary(left, 1e-10));
        EXPECT_TRUE(details::is_unitary(right, 1e-10));
        EXPECT_TRUE(
            details::is_zero(left * Matrix::diagonal(sigma, A.height(), A.width()) * conjugate(right) - A, 1e-10));
    }
}

TEST(SVDTest, TestForLongDoubleMaxSize) {
    using Matrix = Matrix<long double>;

    int operations = 500;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(max_sz, max_sz, max_number);

        Matrix left, right;
        auto sigma = compute_svd(A, &left, &right);

        // check in sigma values > 0
        for (auto s : sigma) {
            EXPECT_TRUE(s >= 0);
        }

        int mx = std::max(A.height(), A.width());

        EXPECT_TRUE(details::is_unitary(left, 1e-10));
        EXPECT_TRUE(details::is_unitary(right, 1e-10));
        EXPECT_TRUE(
            details::is_zero(left * Matrix::diagonal(sigma, A.height(), A.width()) * conjugate(right) - A, 1e-10));
    }
}

TEST(SVDTest, TestForComplex) {
    using Matrix = Matrix<Complex>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(1, max_sz, max_number);

        Matrix left, right;
        auto sigma = compute_svd(A, &left, &right);

        // check in sigma values > 0
        for (auto s : sigma) {
            EXPECT_TRUE(s >= 0);
        }

        std::vector<Complex> sigma_complex;
        for (auto s : sigma) {
            sigma_complex.push_back(Complex(s));
        }

        int mx = std::max(A.height(), A.width());

        EXPECT_TRUE(details::is_unitary(left, 1e-10));
        EXPECT_TRUE(details::is_unitary(right, 1e-10));
        EXPECT_TRUE(details::is_zero(
            left * Matrix::diagonal(sigma_complex, A.height(), A.width()) * conjugate(right) - A, 1e-10));
    }
}

TEST(SVDTest, TestForComplexMaxSize) {
    using Matrix = Matrix<Complex>;

    int operations = 500;
    int max_sz = 100;
    long double max_number = 1e5;

    while (operations-- > 0) {
        auto A = get_random_complex_matrix(max_sz, max_sz, max_number);

        Matrix left, right;
        auto sigma = compute_svd(A, &left, &right);

        // check in sigma values > 0
        for (auto s : sigma) {
            EXPECT_TRUE(s >= 0);
        }

        std::vector<Complex> sigma_complex;
        for (auto s : sigma) {
            sigma_complex.push_back(Complex(s));
        }

        int mx = std::max(A.height(), A.width());

        EXPECT_TRUE(details::is_unitary(left, 1e-10));
        EXPECT_TRUE(details::is_unitary(right, 1e-10));
        EXPECT_TRUE(details::is_zero(
            left * Matrix::diagonal(sigma_complex, A.height(), A.width()) * conjugate(right) - A, 1e-10));
    }
}

TEST(SVDTest, CompareToEigen) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-7;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(1.0, max_sz, max_number);
        Eigen::MatrixXd B(A.height(), A.width());

        for (size_t i = 0; i < A.height(); ++i) {
            for (size_t j = 0; j < A.width(); ++j) {
                B(i, j) = A(i, j);
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(B);
        Vector<long double> eigen_ans(std::min(A.height(), A.width()));

        int ind = 0;
        for (auto elem : svd.singularValues()) {
            eigen_ans[ind++] = elem;
        }

        auto ans = svd_computation::compute_svd(A);

        for (size_t i = 0; i < ans.size(); ++i) {
            EXPECT_NEAR(ans[i], eigen_ans[i], eps);
        }
    }
}

TEST(SVDTest, CompareToEigenMaxSize) {
    using Matrix = Matrix<long double>;

    int operations = 500;
    int max_sz = 100;
    long double max_number = 1e5;

    long double eps = 1e-7;

    while (operations-- > 0) {
        auto A = get_random_double_matrix(max_sz, max_sz, max_number);
        Eigen::MatrixXd B(A.height(), A.width());

        for (size_t i = 0; i < A.height(); ++i) {
            for (size_t j = 0; j < A.width(); ++j) {
                B(i, j) = A(i, j);
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(B);
        Vector<long double> eigen_ans(std::min(A.height(), A.width()));

        int ind = 0;
        for (auto elem : svd.singularValues()) {
            eigen_ans[ind++] = elem;
        }

        auto ans = svd_computation::compute_svd(A);

        for (size_t i = 0; i < ans.size(); ++i) {
            EXPECT_NEAR(ans[i], eigen_ans[i], eps);
        }
    }
}