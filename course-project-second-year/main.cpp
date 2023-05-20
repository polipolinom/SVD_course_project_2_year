#include <iomanip>
#include <iostream>

#include "algorithms/QR_algorithm.h"
#include "algorithms/QR_decomposition.h"
#include "algorithms/bidiagonalization.h"
#include "algorithms/orthonormalize.h"
#include "algorithms/svd_computation.h"
#include "types/matrix.h"
#include "types/vector.h"

int main() {
    using namespace svd_computation;
    Matrix<long double> A = {{1, 2, 3}, {4, 1, 2}, {5, 6, 1}};
    Matrix<long double> V, U;

    auto [Q, R] = get_QR_decomposition(A);
    std::cout << std::setprecision(10) << Q * R << "\n============================\n";
    std::cout << std::setprecision(10) << R << "\n============================\n";
}