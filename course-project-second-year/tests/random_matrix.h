#pragma once
#include "../types/complex.h"
#include "../types/matrix.h"

namespace svd_computation {
Matrix<long double> get_random_double_matrix(int, int, long double);
Matrix<Complex> get_random_complex_matrix(int, int, long double);
}  // namespace svd_computation
