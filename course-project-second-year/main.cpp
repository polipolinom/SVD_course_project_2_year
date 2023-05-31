#include <iomanip>
#include <iostream>
#include <random>

#include "algorithms/QR_algorithm.h"
#include "algorithms/QR_decomposition.h"
#include "algorithms/bidiagonalization.h"
#include "algorithms/orthonormalize.h"
#include "algorithms/svd_computation.h"
#include "types/matrix.h"
#include "types/vector.h"

int main() {
    using namespace svd_computation;
    Matrix<long double> U, V;
    Matrix<long double> A = {{-77509, -20524.6, -81415.3}};

    auto B = compute_svd(A, &V, &U);

    std::cout << U << std::endl;
    Matrix<long double> Sigma(A.height(), A.width());
    for (size_t ind = 0; ind < B.size(); ++ind) {
        Sigma(ind, ind) = B[ind];
    }
    std::cout << std::fixed << std::setprecision(10) << Sigma << "\n\n";
    std::cout << std::setprecision(10) << V * Sigma * transpose(U);

    /*std::default_random_engine gen;
    std::uniform_real_distribution<long double> distribution(-1e5, 1e5);

    while (true) {
        int N = 100;

        Matrix<long double> A(N, N);

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                A(i, j) = distribution(gen);
            }
        }

        auto s = compute_svd(A, &U, &V);
        Matrix<long double> Sigma(N, N);

        for (size_t i = 0; i < N; ++i) {
            Sigma(i, i) = s[i];
        }

        auto ans = U * Sigma * transpose(V) - A;

        long double mx = 0.0;
        long double Mx = 0.0;

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                long double diff = ans(i, j);
                mx = std::max(mx, abs(diff));
                Mx = std::max(Mx, abs(A(i, j)));
            }
        }

        std::cout << mx << " " << Mx << std::endl;
        std::cout << "Done!" << std::endl;
    }*/
}
