#pragma once
#include <utility>

#include "../types/matrix.h"

namespace svd_computation {
inline std::pair<long double, long double> get_givens_rotation(long double a, long double b,
                                                               long double eps = constants::DEFAULT_EPSILON) {
    if (abs(a) < eps) {
        return {1.0, 0.0};
    }
    if (abs(b) > abs(a)) {
        long double t = -a / b;
        return {t / sqrtl(1 + t * t), 1 / sqrtl(1 + t * t)};
    }
    long double t = -b / a;
    return {1 / sqrtl(1 + t * t), t / sqrtl(1 + t * t)};
}

template <typename Type>
inline void multiply_right_givens(Matrix<Type>& A, Type c, Type s, int i, int j) {
    for (size_t k = 0; k < A.height(); ++k) {
        Type a = A(k, i);
        Type b = A(k, j);

        A(k, i) = c * a - s * b;
        A(k, j) = s * a + c * b;
    }
}

template <typename Type>
inline void multiply_left_givens(Matrix<Type>& A, Type c, Type s, int i, int j) {
    for (size_t k = 0; k < A.width(); ++k) {
        Type a = A(i, k);
        Type b = A(j, k);

        A(i, k) = c * a - s * b;
        A(j, k) = s * a + c * b;
    }
}
}  // namespace svd_computation
