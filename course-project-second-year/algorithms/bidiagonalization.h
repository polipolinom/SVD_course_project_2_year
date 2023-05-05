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
void left_reflection(Matrix<Type>& A, size_t ind, Matrix<Type>* left_basis = nullptr,
                     long double eps = constants::DEFAULT_EPSILON) {
    set_low_values_zero(A, eps);
    Matrix<Type> u(A.height(), 1);
    long double s = 0.0;
    for (size_t k = ind; k < A.height(); k++) {
        s += abs(A(k, ind)) * abs(A(k, ind));
    }
    s = sqrtl(s);

    // all numbers less than eps sets to zero with function set_low_values_zero
    if (s == 0.0) {
        return;
    }

    Type alpha;
    if (A(ind, ind) == 0.0) {
        alpha = Type(s);
    } else {
        alpha = s * A(ind, ind) / abs(A(ind, ind));
    }

    u(ind, 0) = sqrt((Type(1.0) + abs(A(ind, ind)) / s) / Type(2.0));
    Type coef = Type(1.0) / (Type(2.0) * alpha * u(ind, 0));
    for (size_t k = ind + 1; k < A.height(); k++) {
        u(k, 0) = coef * A(k, ind);
    }

    A -= Type(2.0) * u * (u.conjugate() * A);
    if (left_basis != nullptr) {
        (*left_basis) -= Type(2.0) * u * u.conjugate() * (*left_basis);
    }
}

template <typename Type>
void right_reflection(Matrix<Type>& A, size_t ind, Matrix<Type>* right_basis = nullptr,
                      long double eps = constants::DEFAULT_EPSILON) {
    set_low_values_zero(A, eps);
    Matrix<Type> u(A.width(), 1);
    long double s = 0;
    for (size_t k = ind + 1; k < A.width(); k++) {
        s += abs(A(ind, k)) * abs(A(ind, k));
    }
    s = sqrtl(s);

    // all numbers less than eps sets to zero with function set_low_values_zero
    if (s == 0.0) {
        return;
    }

    Type alpha;
    if (A(ind, ind + 1) == 0.0) {
        alpha = Type(s);
    } else {
        alpha = s * A(ind, ind + 1) / abs(A(ind, ind + 1));
    }

    u(ind + 1, 0) = sqrt((Type(1.0) + abs(A(ind, ind + 1)) / s) / Type(2.0));
    Type coef = Type(1.0) / (Type(2.0) * alpha * u(ind + 1, 0));
    for (size_t k = ind + 2; k < A.width(); ++k) {
        u(k, 0) = coef * conjugate(A(ind, k));
    }

    A -= Type(2.0) * (A * u) * u.conjugate();
    if (right_basis != nullptr) {
        (*right_basis) -= Type(2.0) * (*right_basis) * u * u.conjugate();
    }
}
}  // namespace details

template <typename Type>
Matrix<Type> bidiagonalize(Matrix<Type>& A, Matrix<Type>* left_basis = nullptr, Matrix<Type>* right_basis = nullptr,
                           long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> B = A;
    if (left_basis != nullptr) {
        *left_basis = Matrix<Type>::ones(A.height());
    }
    if (right_basis != nullptr) {
        *right_basis = Matrix<Type>::ones(A.width());
    }
    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        details::left_reflection(B, ind, left_basis, eps);
        if (ind + 1 < A.width()) {
            details::right_reflection(B, ind, right_basis, eps);
        }
    }
    set_low_values_zero(B, eps);
    return B;
}
}  // namespace svd_computation