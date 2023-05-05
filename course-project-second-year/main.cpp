#include <iostream>

#include "algorithms/QR_decomposition.h"
#include "algorithms/bidiagonalization.h"
#include "algorithms/orthonormalize.h"
#include "types/matrix.h"
#include "types/vector.h"

int main() {
    using namespace svd_computation;
    Matrix<Complex> A = {{2, Complex(3, 1), 4, 5, 6}, {1, 2, Complex(3, 1), 4, 5}};
    Matrix<Complex> U, V;

    auto B = bidiagonalize(A, &U, &V);

    std::cout << B << "\n============================\n";
    std::cout << U.conjugate() * U << "\n============================\n";
    std::cout << V.conjugate() * V << "\n============================\n";
}