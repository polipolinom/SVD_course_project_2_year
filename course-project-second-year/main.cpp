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
    Matrix<long double> A = {{799.408, -519.941, -359.672, 599.85},
                             {498.816, -399.882, -449.345, 747.7},
                             {2498.42, -1599.84, -1049.13, 1749.6},
                             {2498, -1599.88, -1049.34, 1749.7}};
    Matrix<long double> V, U;

    auto B = compute_svd(A, &V, &U);
    std::cout << std::setprecision(10) << V * B * transpose(U) - A << "\n============================\n";
}
