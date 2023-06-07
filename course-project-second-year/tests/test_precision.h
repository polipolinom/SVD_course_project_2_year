#pragma once

#include <gtest/gtest.h>

#include <chrono>
#include <ctime>
#include <random>

#include "../algorithms/svd_computation.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../types/vector.h"
#include "Eigen/eigen"
#include "random_objects.h"

using namespace svd_computation;
using namespace std::chrono;

std::vector<long double> test_precision_eigen(const std::vector<size_t> &ns) {
    long double max_number = 1e5;
    size_t iterations_count = 10;

    std::vector<long double> res;
    for (size_t n : ns) {
        long double total_diff = 0;
        for (size_t i = 0; i < iterations_count; i++) {
            Matrix<long double> A = get_random_double_matrix(n, n, max_number);
            std::vector<long double> my_sigma = svd_computation::compute_svd(A);

            Eigen::MatrixXd B(n, n);
            for (size_t x = 0; x < n; x++) {
                for (size_t y = 0; y < n; y++) {
                    B(x, y) = A(x, y);
                }
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> eigen_sigma(B);
            svd_computation::compute_svd(A);

            long double max_diff = 0;
            size_t ptr = 0;
            for (auto elem : eigen_sigma.singularValues()) {
                max_diff = std::max(max_diff, abs(elem - my_sigma[ptr++]));
            }

            total_diff += max_diff;
        }
        res.push_back(total_diff / iterations_count);
    }
    return res;
}

std::vector<long double> test_precision(const std::vector<size_t> &ns) {
    using Matrix = Matrix<long double>;

    long double max_number = 1e5;
    size_t iterations_count = 10;

    std::vector<long double> res;
    for (size_t n : ns) {
        long double total_diff = 0;
        for (size_t i = 0; i < iterations_count; i++) {
            Matrix A = get_random_double_matrix(n, n, max_number);
            Matrix V, U;
            std::vector sigma = svd_computation::compute_svd(A, &V, &U);

            auto diff = V * Matrix::diagonal(sigma, A.height(), A.width()) * transpose(U) - A;

            long double max_diff = 0;

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    max_diff = std::max(diff(i, j), max_diff);
                }
            }

            total_diff += max_diff;
        }
        res.push_back(total_diff / iterations_count);
    }
    return res;
}

// for work need to delete constexpr у constants::MAX_OPERATIONS and add in CmakeList executable for file
// tests/graphic_tests.cpp
std::vector<long double> test_precision_operations(const std::vector<size_t> &operations) {
    using Matrix = Matrix<long double>;

    long double max_number = 1e5;
    size_t iterations_count = 100;
    size_t n = 100;

    std::vector<long double> res;
    for (size_t operation : operations) {
        constants::MAX_OPERATIONS = operation;

        long double total_diff = 0;
        for (size_t i = 0; i < iterations_count; i++) {
            Matrix A = get_random_double_matrix(n, n, max_number);
            Matrix V, U;
            std::vector sigma = svd_computation::compute_svd(A, &V, &U);

            auto diff = V * Matrix::diagonal(sigma, A.height(), A.width()) * transpose(U) - A;

            long double max_diff = 0;

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    max_diff = std::max(diff(i, j), max_diff);
                }
            }

            total_diff += max_diff;
        }
        res.push_back(total_diff / iterations_count);
    }

    constants::MAX_OPERATIONS = 50;
    return res;
}

// for work need to delete constexpr у constants::DEFAULT_EPSILON and add in CmakeList executable for file
// tests/graphic_tests.cpp
std::vector<long double> test_precision_epsilons(const std::vector<long double> &epsilons) {
    using Matrix = Matrix<long double>;

    long double max_number = 1e5;
    size_t iterations_count = 100;
    size_t n = 100;

    std::vector<long double> res;
    for (long double eps : epsilons) {
        constants::DEFAULT_EPSILON = eps;

        long double total_diff = 0;
        for (size_t i = 0; i < iterations_count; i++) {
            Matrix A = get_random_double_matrix(n, n, max_number);
            Matrix V, U;
            std::vector sigma = svd_computation::compute_svd(A, &V, &U);

            auto diff = V * Matrix::diagonal(sigma, A.height(), A.width()) * transpose(U) - A;

            long double max_diff = 0;

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    max_diff = std::max(diff(i, j), max_diff);
                }
            }

            total_diff += max_diff;
        }
        res.push_back(total_diff / iterations_count);
    }

    constants::DEFAULT_EPSILON = 1e-16;
    return res;
}
