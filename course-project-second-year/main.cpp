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
    Matrix<long double> A = {{1, 2, 3}, {0, 1, 2}, {0, 0, 1}};
    Matrix<long double> V, U;

    auto B = compute_svd(A, &V, &U);
    std::cout << std::setprecision(10) << B << "\n============================\n";
}