#pragma once

#include "../algorithms/constants.h"
#include "../types/matrix.h"

namespace svd_computation {
namespace details {
template <typename Type>
void set_low_values_zero(Matrix<Type>& A, const long double eps = constants::DEFAULT_EPSILON) {
    for (size_t i = 0; i < A.height(); ++i) {
        for (size_t j = 0; j < A.width(); ++j) {
            if (abs(A(i, j)) <= eps) {
                A(i, j) = Type(0.0);
            }
        }
    }
}
}  // namespace details
}  // namespace svd_computation