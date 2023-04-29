#include <iostream>

#include "algorithms/QR_decomposition.h"
#include "algorithms/orthonormalize.h"
#include "types/matrix.h"
#include "types/vector.h"

int main() {
    using namespace svd_computation;
    Matrix<long double> A = {{2, 3}, {2, 3}};

    auto [Q, R] = get_QR_decomposition(A);

    std::cout << Q << std::endl << "========================================" << std::endl;
    std::cout << R << std::endl << "========================================" << std::endl;
    std::cout << Q.transpose() * Q << std::endl << "========================================" << std::endl;
    std::cout << Q * R << std::endl << "========================================" << std::endl;
}