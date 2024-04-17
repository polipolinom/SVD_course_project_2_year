#include <iostream>

#include "algorithms/svd_computation.h"
#include "types/matrix.h"

using namespace svd_computation;

int main() {
    Matrix<long double> A = {
        {0.177341, -2.55219, -0.00371925}, {0.177341, -2.55219, -0.00371925}, {0, 0.00265251, -2.56133}};
    auto diag = compute_svd(A);
    for (auto i : diag) {
        std::cout << i << " ";
    }
    std::cout << (A == A) << std::endl;
}
