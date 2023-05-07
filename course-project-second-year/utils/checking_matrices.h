#pragma once

#include "../algorithms/constants.h"
#include "../types/matrix.h"

namespace svd_computation {
namespace details {
template <typename Type>
bool is_upper_triangular(const Matrix<Type>& A, const long double eps = constants::DEFAULT_EPSILON) {
    for (size_t i = 0; i < A.height(); ++i) {
        for (size_t j = 0; j < std::min(i, A.width()); ++j) {
            if (abs(A(i, j)) > eps) {
                return false;
            }
        }
    }
    return true;
}

template <typename Type>
bool is_tridiagonal(const Matrix<Type>& A, const long double eps = constants::DEFAULT_EPSILON) {
    if (A.height() != A.width()) {
        return false;
    }
    for (size_t i = 0; i < A.height(); ++i) {
        for (size_t j = 0; j + 1 < i; ++j) {
            if (abs(A(i, j)) > eps) {
                return false;
            }
        }
        for (size_t j = i + 2; j < A.width(); ++j) {
            if (abs(A(i, j)) > eps) {
                return false;
            }
        }
    }
    return true;
}

template <typename Type>
bool is_diagonal(const Matrix<Type>& A, const long double eps = constants::DEFAULT_EPSILON) {
    if (A.height() != A.width()) {
        return false;
    }
    for (size_t i = 0; i < A.height(); ++i) {
        for (size_t j = 0; j < A.width(); ++j) {
            if (i != j && abs(A(i, j)) > eps) {
                return false;
            }
        }
    }
    return true;
}
}  // namespace details
}  // namespace svd_computation