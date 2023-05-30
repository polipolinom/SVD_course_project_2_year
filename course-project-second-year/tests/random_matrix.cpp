#include "random_matrix.h"

#include <cassert>
#include <chrono>
#include <ctime>
#include <random>

namespace svd_computation {
Matrix<long double> get_random_double_matrix(int min_sz, int max_sz, long double max_number) {
    assert(max_sz > 0);
    assert(min_sz > 0);

    std::mt19937 rnd(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<long double> distribution(-max_number, max_number);

    int n = rnd() % (max_sz - min_sz + 1) + min_sz;
    int m = rnd() % (max_sz - min_sz + 1) + min_sz;

    Matrix<long double> A(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A(i, j) = distribution(rnd);
        }
    }

    return A;
}

Matrix<Complex> get_random_complex_matrix(int min_sz, int max_sz, long double max_number) {
    assert(max_sz > 0);
    assert(min_sz > 0);

    std::mt19937 rnd(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<long double> distribution(-max_number, max_number);

    int n = rnd() % (max_sz - min_sz + 1) + min_sz;
    int m = rnd() % (max_sz - min_sz + 1) + min_sz;

    Matrix<Complex> A(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A(i, j) = Complex(distribution(rnd), distribution(rnd));
        }
    }

    return A;
}
}  // namespace svd_computation