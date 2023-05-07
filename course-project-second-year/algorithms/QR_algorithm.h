#pragma once

#include <cassert>

#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "../utils/matrix_split_join.h"
#include "../utils/set_values.h"
#include "QR_decomposition.h"
#include "constants.h"

namespace svd_computation {
template <typename Type>
Matrix<Type> get_Schur_decomposition(const Matrix<Type>& A, const Type shift = 0.0, Matrix<Type>* basis = nullptr,
                                     const long double eps = constants::DEFAULT_EPSILON) {
    assert(A.height() == A.width());

    Matrix<Type> result = A;
    if (basis != nullptr) {
        (*basis) = Matrix<Type>::ones(A.height());
    }

    while (!details::is_upper_triangular(result, eps)) {
        details::set_low_values_zero(result, eps);
        auto result_with_shift = result - shift * Matrix<Type>::ones(result.height());
        auto [Q, R] = get_QR_decomposition(result_with_shift, eps);
        result = R * Q + shift * Matrix<Type>::ones(result.height());
        if (basis != nullptr) {
            (*basis) *= Q;
        }
    }
    return result;
}

namespace details {
long double get_Wilkinson_shift(const Matrix<long double>& A) {
    assert(A.height() > 1);
    size_t n = A.height();
    long double shift = A(n - 1, n - 1);
    long double sigma = (A(n - 2, n - 2) - A(n - 1, n - 1)) / 2.0;

    long double coef = -1.0;

    if (sigma < 0) {
        coef = 1.0;
    }
    shift += coef * (A(n - 1, n - 2) * A(n - 1, n - 2)) /
             (abs(sigma) + sqrtl(sigma * sigma + (A(n - 1, n - 2) * A(n - 1, n - 2))));

    return shift;
}

Matrix<long double> apply_qr_for_tridiagonal(const Matrix<long double>& A, Matrix<long double>* basis,
                                             const long double eps = constants::DEFAULT_EPSILON) {
    assert(is_tridiagonal(A, eps));
    if (A.height() == 1) {
        if (basis != nullptr) {
            *basis = {{1.0}};
        }
        return A;
    }
    Matrix<long double> result = A;
    while (!is_diagonal(result, eps)) {
        set_low_values_zero(result, eps);
        for (size_t ind = 0; ind + 1 < A.height(); ++ind) {
            if (A(ind, ind + 1) == 0.0) {
                auto [result1, result2] = split_matrix(result, ind, ind);
                Matrix<long double> basis1 = Matrix<long double>::ones(ind + 1);
                Matrix<long double> basis2 = Matrix<long double>::ones(A.height() - ind - 1);
                result1 = apply_qr_for_tridiagonal(result1, &basis1, eps);
                result2 = apply_qr_for_tridiagonal(result2, &basis2, eps);
                result = join_matrix(result1, result2);
                if (basis != nullptr) {
                    (*basis) *= join_matrix(basis1, basis2);
                }
                return result;
            }
        }
        long double shift = get_Wilkinson_shift(result);
        auto result_with_shift = result - shift * Matrix<long double>::ones(result.height());
        auto [Q, R] = get_QR_decomposition(result_with_shift, eps);
        result = R * Q + shift * Matrix<long double>::ones(result.height());
        if (basis != nullptr) {
            (*basis) *= Q;
        }
    }
    set_low_values_zero(result, eps);
    return result;
}
}  // namespace details
}  // namespace svd_computation