#pragma once

#include "../types/matrix.h"

namespace svd_computation {
namespace details {
template <typename Type>
long double column_abs_under(const Matrix<Type>& A, const size_t column, const size_t first_ind) {
    assert(column >= 0 && column < A.width());
    long double s = 0.0;
    for (size_t k = first_ind; k < A.height(); ++k) {
        s += abs(A(k, column)) * abs(A(k, column));
    }
    s = sqrtl(s);
    return s;
}

template <typename Type>
long double row_abs_under(const Matrix<Type>& A, const size_t row, const size_t first_ind) {
    assert(row >= 0 && row < A.height());
    long double s = 0.0;
    for (size_t k = first_ind; k < A.width(); ++k) {
        s += abs(A(row, k)) * abs(A(row, k));
    }
    s = sqrtl(s);
    return s;
}
}  // namespace details

template <typename Type>
long double left_reflection(Matrix<Type>& A, const int row, const int column, Matrix<Type>* left_basis = nullptr,
                            const long double eps = constants::DEFAULT_EPSILON) {
    assert(row >= 0 && row < A.height());
    assert(column >= 0 && column < A.width());

    // details::set_low_values_zero(A, eps);
    Matrix<Type> u(A.height(), 1);

    long double s = details::column_abs_under(A, row, column);

    if (s <= eps) {
        return 0.0;
    }

    Type alpha = Type(s);
    if (abs(A(row, column)) > eps) {
        alpha *= A(row, column) / abs(A(row, column));
    }

    u(row, 0) = A(row, column) + alpha;
    for (size_t k = row + 1; k < A.height(); ++k) {
        u(k, 0) = A(k, column);
    }

    long double coef = details::column_abs_under(u, 0, row);
    // coef > 0 because exist |A(k, ind)| > eps / A.height() (otherwise s <= eps)
    u /= coef;

    Matrix<Type> P = Matrix<Type>::identity(A.height()) - Type(2.0) * u * conjugate(u);

    A = P * A;
    if (left_basis != nullptr) {
        (*left_basis) = P * (*left_basis);
    }

    return -alpha;
}

template <typename Type>
long double right_reflection(Matrix<Type>& A, const int row, const int column, Matrix<Type>* right_basis = nullptr,
                             const long double eps = constants::DEFAULT_EPSILON) {
    assert(row >= 0 && row < A.height());
    assert(column >= 0 && column < A.width());

    // details::set_low_values_zero(A, eps);
    Matrix<Type> u(1, A.width());
    long double s = details::row_abs_under(A, row, column);

    if (s <= eps) {
        return 0.0;
    }

    Type alpha = Type(s);
    if (abs(A(row, column)) > eps) {
        alpha *= A(row, column) / abs(A(row, column));
    }

    u(0, column) = A(row, column) + alpha;
    for (size_t k = column + 1; k < A.width(); ++k) {
        u(0, k) = A(row, k);
    }

    long double coef = details::row_abs_under(u, 0, column);
    // coef > 0 because exist |A(ind, k)| > eps / A.width() (otherwise s <= eps)
    u /= coef;

    Matrix<Type> P = Matrix<Type>::identity(A.width()) - Type(2.0) * conjugate(u) * u;

    A *= P;

    if (right_basis != nullptr) {
        (*right_basis) *= P;
    }

    return -alpha;
}
}  // namespace svd_computation
