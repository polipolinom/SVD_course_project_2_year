#pragma once

#include <cmath>
#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/complement_orthobase.h"
#include "../utils/set_values.h"

namespace svd_computation {
namespace details {
template <typename Type>
long double column_abs(const Matrix<Type>& A, const size_t column, const size_t first_ind) {
    assert(column >= 0 && column < A.width());
    long double s = 0.0;
    for (size_t k = first_ind; k < A.height(); ++k) {
        s += abs(A(k, column)) * abs(A(k, column));
    }
    s = sqrtl(s);
    return s;
}

template <typename Type>
long double row_abs(const Matrix<Type>& A, const size_t row, const size_t first_ind) {
    assert(row >= 0 && row < A.height());
    long double s = 0.0;
    for (size_t k = first_ind; k < A.width(); ++k) {
        s += abs(A(row, k)) * abs(A(row, k));
    }
    s = sqrtl(s);
    return s;
}

template <typename Type>
long double left_reflection(Matrix<Type>& A, const size_t ind, Matrix<Type>* left_basis = nullptr,
                            const long double eps = constants::DEFAULT_EPSILON) {
    assert(ind >= 0 && ind < A.height());
    details::set_low_values_zero(A, eps);
    Matrix<Type> u(A.height(), 1);

    long double s = column_abs(A, ind, ind);

    // all numbers less than eps sets to zero with function set_low_values_zero
    if (s == 0.0) {
        return 0.0;
    }

    Type alpha = Type(s);
    if (A(ind, ind) != 0.0) {
        alpha *= A(ind, ind) / abs(A(ind, ind));
    }

    u(ind, 0) = A(ind, ind) + alpha;
    for (size_t k = ind + 1; k < A.height(); ++k) {
        u(k, 0) = A(k, ind);
    }

    long double coef = column_abs(u, 0, ind);
    u /= coef;

    Matrix<Type> P = Matrix<Type>::ones(A.height()) - Type(2.0) * u * u.conjugate();

    A = P * A;
    if (left_basis != nullptr) {
        (*left_basis) = P * (*left_basis);
    }

    return -alpha;
}

template <typename Type>
long double right_reflection(Matrix<Type>& A, const size_t ind, Matrix<Type>* right_basis = nullptr,
                             const long double eps = constants::DEFAULT_EPSILON) {
    assert(ind >= 0 && ind + 1 < A.width());
    details::set_low_values_zero(A, eps);
    Matrix<Type> u(A.width(), 1);
    long double s = row_abs(A, ind, ind + 1);

    // all numbers less than eps sets to zero with function set_low_values_zero
    if (s == 0.0) {
        return 0.0;
    }

    Type alpha = Type(s);
    if (A(ind, ind + 1) != 0.0) {
        alpha *= A(ind, ind + 1) / abs(A(ind, ind + 1));
    }

    u(ind + 1, 0) = A(ind, ind + 1) + alpha;
    for (size_t k = ind + 2; k < A.width(); ++k) {
        u(k, 0) = A(ind, k);
    }

    long double coef = column_abs(u, 0, ind + 1);
    u /= coef;

    Matrix<Type> P = Matrix<Type>::ones(A.width()) - Type(2.0) * u.conjugate().transpose() * u.transpose();

    A *= P;

    if (right_basis != nullptr) {
        (*right_basis) *= P;
    }

    return -alpha;
}
}  // namespace details

template <typename Type>
Matrix<long double> bidiagonalize(const Matrix<Type>& A, Matrix<Type>* left_basis = nullptr,
                                  Matrix<Type>* right_basis = nullptr,
                                  const long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> B = A;
    Matrix<long double> result(A.height(), A.width());
    if (left_basis != nullptr) {
        *left_basis = Matrix<Type>::ones(A.height());
    }
    if (right_basis != nullptr) {
        *right_basis = Matrix<Type>::ones(A.width());
    }
    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        result(ind, ind) = details::left_reflection(B, ind, left_basis, eps);
        if (ind + 1 < A.width()) {
            result(ind, ind + 1) = details::right_reflection(B, ind, right_basis, eps);
        }
    }
    if (left_basis != nullptr) {
        (*left_basis) = (*left_basis).conjugate();
    }
    return result;
}
}  // namespace svd_computation