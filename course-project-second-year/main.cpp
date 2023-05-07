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
    Matrix<long double> A = {{2, 3}, {4, 5}, {5, 6}, {7, 6}};
    Matrix<long double> U, V;

    auto B = compute_svd(A, &U, &V);

    std::cout << B << "\n\n";
    std::cout << U * B * V.transpose();
}