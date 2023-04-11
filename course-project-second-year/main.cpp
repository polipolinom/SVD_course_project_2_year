#include <iostream>

#include "types/complex.h"
#include "types/matrix.h"

using namespace svd_computation;

signed main() {
    Matrix<Complex> A(Complex(2), 5);
    Matrix<Complex> B(Complex(3), 5);
    Matrix<Complex> C = Complex(1, 2) * A;
    std::cout << A << std::endl << B << std::endl << C << std::endl;
    return 0;
}