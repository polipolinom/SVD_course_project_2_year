#pragma once

#include "../types/complex.h"
#include "../types/matrix.h"
#include "../types/vector.h"

namespace svd_computation {
Matrix<long double> get_random_double_matrix(int, int, long double);
Matrix<Complex> get_random_complex_matrix(int, int, long double);
Vector<long double> get_random_double_vector(int, int, long double);
Vector<Complex> get_random_complex_vector(int, int, long double);
}  // namespace svd_computation
