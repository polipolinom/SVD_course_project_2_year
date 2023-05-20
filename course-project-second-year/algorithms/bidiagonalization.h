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
long double left_reflection(Matrix<Type>& A, const size_t ind, Matrix<Type>* left_basis = nullptr,
                            const long double eps = constants::DEFAULT_EPSILON) {
    assert(ind >= 0 && ind < A.height());
    details::set_low_values_zero(A, eps);
    Matrix<Type> u(A.height(), 1);

    long double s = column_abs_under(A, ind, ind);

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

    long double coef = column_abs_under(u, 0, ind);
    // coef > 0 because exist |A(k, ind)| > eps (otherwise s = 0.0)
    u /= coef;

    Matrix<Type> P = Matrix<Type>::identity(A.height()) - Type(2.0) * u * conjugate(u);

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
    long double s = row_abs_under(A, ind, ind + 1);

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

    long double coef = column_abs_under(u, 0, ind + 1);
    // coef > 0 because exist |A(ind, k)| > eps (otherwise s = 0.0)
    u /= coef;

    u.transpose();

    Matrix<Type> P = Matrix<Type>::identity(A.width()) - Type(2.0) * conjugate(u) * u;

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
        *left_basis = Matrix<Type>::identity(A.height());
    }
    if (right_basis != nullptr) {
        *right_basis = Matrix<Type>::identity(A.width());
    }
    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        result(ind, ind) = details::left_reflection(B, ind, left_basis, eps);
        if (ind + 1 < A.width()) {
            result(ind, ind + 1) = details::right_reflection(B, ind, right_basis, eps);
        }
    }
    if (left_basis != nullptr) {
        (*left_basis).conjugate();
    }
    return result;
}

namespace details {
Matrix<long double> bidiagonalize_with_right_basis(Matrix<long double>& A, const Matrix<long double>& right_basis,
                                                   const long double eps) {
    using Matrix = Matrix<long double>;

    Matrix basis = Matrix::identity(A.height());

    Matrix A1 = A * right_basis;

    long double alpha = 0;
    long double beta = 0;
    Vector<long double> u(0, A.height());
    for (size_t ind = 0; ind < A.height(); ++ind) {
        alpha = abs(A * right_basis.column(ind) - beta * u);
        u = (A1.column(ind) - beta * u) / alpha;
        for (size_t i = 0; i < u.size(); ++i) {
            basis(i, ind) = u[i];
        }

        if (ind + 1 < A.width()) {
            beta = dot_product(transpose(u), A1.column(ind + 1));
        }
    }
    std::cout << "#######################\n";
    std::cout << transpose(basis) * A * right_basis << "\n=========================\n";
    return basis;
}
}  // namespace details
}  // namespace svd_computation