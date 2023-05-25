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
    Matrix<long double> A = {{-999.998, 200., -899.999, 1500.},
                             {-4000., 1400., -1800., 3000.},
                             {-2000., 200., -2400., 4000.},
                             {-2000., 200., -2400., 4000.}};
    Matrix<long double> U, V;

    // A = {{1, 1, 0}, {0, -1, 1}, {0, 0, 1}};
    auto B = compute_svd(A, &V, &U);
    std::cout << std::setprecision(10) << V * B * transpose(U);

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
