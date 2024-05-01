#include <iostream>

#include "algorithms/svd_computation.h"
#include "types/matrix.h"

using namespace svd_computation;

int main() {
    Matrix<long double> A = {{0.177341, -2.55219, 0}, {0, -2.55219, -0.00371925}, {0, 0, -2.56133}};
    std::cout << transpose(A) * A << std::endl << std::endl;
    std::cout << mult_band(transpose(A), 2, 1, A, 1, 2) << std::endl;
}
