#pragma once

#include "../types/complex.h"
#include "../types/matrix.h"
#include "constants.h"

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

template <typename Type>
void set_value_real_left(Matrix<Type>& A, const int row, const int column, Matrix<Type>* left_basis = nullptr,
                         const long double eps = constants::DEFAULT_EPSILON) {
    if (abs(A(row, column)) <= eps) {
        return;
    }

    Type coef = conjugate(A(row, column) / abs(A(row, column)));
    for (size_t i = 0; i < A.width(); i++) {
        A(row, i) *= coef;
    }

    if (left_basis != nullptr) {
        for (size_t i = 0; i < A.height(); ++i) {
            (*left_basis)(row, i) *= coef;
        }
    }
}

template <typename Type>
void set_value_real_right(Matrix<Type>& A, const int row, const int column, Matrix<Type>* right_basis = nullptr,
                          const long double eps = constants::DEFAULT_EPSILON) {
    if (abs(A(row, column)) <= eps) {
        return;
    }
    Type coef = conjugate(A(row, column) / abs(A(row, column)));
    for (size_t i = 0; i < A.height(); i++) {
        A(i, column) *= coef;
    }
    if (right_basis != nullptr) {
        for (size_t i = 0; i < A.width(); ++i) {
            (*right_basis)(i, column) *= coef;
        }
    }
}
}  // namespace details

template <typename Type>
void left_reflection(Matrix<Type>& A, const int row, const int column, Matrix<Type>* left_basis = nullptr,
                     const long double eps = constants::DEFAULT_EPSILON) {
    assert(row >= 0 && row < A.height());
    assert(column >= 0 && column < A.width());

    Matrix<Type> u(A.height(), 1);

    long double s = details::column_abs_under(A, column, row);

    if (s <= eps) {
        details::set_value_real_left(A, row, column, left_basis, eps);
        return;
    }

    Type alpha = Type(s);
    if (abs(A(row, column)) > eps) {
        alpha *= A(row, column) / abs(A(row, column));
    }

    u(row, 0) = A(row, column) - alpha;
    for (size_t k = row + 1; k < A.height(); ++k) {
        u(k, 0) = A(k, column);
    }

    long double coef = details::column_abs_under(u, 0, 0);
    if (coef <= eps) {
        details::set_value_real_left(A, row, column, left_basis, eps);
        return;
    }
    u /= coef;

    A -= Type(2.0) * u * (conjugate(u) * A);
    if (left_basis != nullptr) {
        (*left_basis) -= Type(2.0) * u * (conjugate(u) * (*left_basis));
    }

    details::set_value_real_left(A, row, column, left_basis, eps);
}

template <typename Type>
void right_reflection(Matrix<Type>& A, const int row, const int column, Matrix<Type>* right_basis = nullptr,
                      const long double eps = constants::DEFAULT_EPSILON) {
    assert(row >= 0 && row < A.height());
    assert(column >= 0 && column < A.width());

    Matrix<Type> u(1, A.width());
    long double s = details::row_abs_under(A, row, column);

    if (s <= eps) {
        details::set_value_real_right(A, row, column, right_basis, eps);
        return;
    }

    Type alpha = Type(s);
    if (abs(A(row, column)) > eps) {
        alpha *= A(row, column) / abs(A(row, column));
    }

    u(0, column) = A(row, column) - alpha;
    for (size_t k = column + 1; k < A.width(); ++k) {
        u(0, k) = A(row, k);
    }

    long double coef = details::row_abs_under(u, 0, 0);
    if (coef <= eps) {
        details::set_value_real_right(A, row, column, right_basis, eps);
        return;
    }

    u /= coef;

    A -= Type(2.0) * (A * conjugate(u)) * u;

    if (right_basis != nullptr) {
        (*right_basis) -= Type(2.0) * ((*right_basis) * conjugate(u)) * u;
    }

    details::set_value_real_right(A, row, column, right_basis, eps);
    return;
}
}  // namespace svd_computation
