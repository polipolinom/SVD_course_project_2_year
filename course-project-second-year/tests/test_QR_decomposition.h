#include <gtest/gtest.h>

#include <random>

#include "../algorithms/QR_decomposition.h"
#include "../types/complex.h"
#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/checking_matrices.h"

using namespace svd_computation;

TEST(QRDecompositionTest, QRForLongDouble) {
    using Matrix = Matrix<long double>;

    int operations = 1000;
    int max_sz = 100;
    long double max_number = 1e5;

    std::default_random_engine gen;
    std::uniform_real_distribution<long double> distribution(-max_number, max_number);

    std::mt19937 random_int;

    while (operations-- > 0) {
        int n = random_int() % max_sz + 1;
        int m = random_int() % max_sz + 1;

        Matrix A(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                A(i, j) = distribution(gen);
            }
        }

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

    std::default_random_engine gen;
    std::uniform_real_distribution<long double> distribution(-max_number, max_number);

    std::mt19937 random_int;

    while (operations-- > 0) {
        int n = random_int() % max_sz + 1;
        int m = random_int() % max_sz + 1;

        Matrix A(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                A(i, j) = Complex(distribution(gen), distribution(gen));
            }
        }

        auto [Q, R] = get_QR_decomposition(A);

        EXPECT_TRUE(details::is_upper_triangular(R));
        EXPECT_TRUE(details::is_unitary(Q));
        EXPECT_TRUE(details::is_zero(Q * R - A));
    }
}
