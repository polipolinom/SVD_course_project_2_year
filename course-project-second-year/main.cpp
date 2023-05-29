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
    Matrix<long double> A = {
        {-999.998, 200., -899.999}, {-4000., 1400., -1800.}, {-2000., 200., -2400.}, {-2000., 200., -2400.}};
    Matrix<long double> U, V;

    // A = {{1, 1, 0}, {0, -1, 1}, {0, 0, 1}};
    auto B = compute_svd(A, &V, &U);
    Matrix<long double> Sigma(A.height(), A.width());
    for (size_t ind = 0; ind < B.size(); ++ind) {
        Sigma(ind, ind) = B[ind];
    }
    std::cout << std::setprecision(10) << V * Sigma * transpose(U);

    /*std::default_random_engine gen;
    std::uniform_real_distribution<long double> distribution(0.0, 1.0);

    while (true) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                A(i, j) = distribution(gen);
            }
        }
        std::cout << compute_svd(A) << "\n";
        std::cout << "Done!" << std::endl;
    }*/
}
