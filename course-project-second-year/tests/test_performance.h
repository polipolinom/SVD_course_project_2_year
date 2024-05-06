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

std::vector<std::pair<int, int>> test_performance(const std::vector<size_t> &ns) {
    long double max_number = 1e5;
    size_t iterations_count = 10;

    std::vector<std::pair<int, int>> res;
    for (size_t n : ns) {
        int total_time = 0;
        int total_time_eigen = 0;
        for (size_t i = 0; i < iterations_count; i++) {
            Matrix<long double> A = get_random_double_matrix(n, n, max_number);
            auto start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            svd_computation::compute_svd(A);
            auto finish = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            total_time += finish - start;

            Eigen::MatrixXd B(n, n);
            for (size_t x = 0; x < n; x++) {
                for (size_t y = 0; y < n; y++) {
                    B(x, y) = A(x, y);
                }
            }
            start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            Eigen::JacobiSVD<Eigen::MatrixXd> eigen_sigma(B);
            finish = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            eigen_sigma.singularValues();
            total_time_eigen += finish - start;
        }
        res.emplace_back(total_time / iterations_count, total_time_eigen / iterations_count);
    }
    return res;
}
