#pragma once

#include <cassert>

#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "../utils/matrix_split_join.h"
#include "../utils/set_values.h"
#include "QR_decomposition.h"
#include "bidiagonalization.h"
#include "constants.h"

namespace svd_computation {
template <typename Type>
Matrix<Type> get_Schur_decomposition(const Matrix<Type>& A, const Type shift = 0.0, Matrix<Type>* basis = nullptr,
                                     const long double eps = constants::DEFAULT_EPSILON) {
    assert(A.height() == A.width());

    Matrix<Type> result = A;
    if (basis != nullptr) {
        (*basis) = Matrix<Type>::identity(A.height());
    }

    while (!details::is_upper_triangular(result, eps)) {
        details::set_low_values_zero(result, eps);
        auto result_with_shift = result - shift * Matrix<Type>::identity(result.height());
        auto [Q, R] = get_QR_decomposition(result_with_shift, eps);
        result = R * Q + shift * Matrix<Type>::identity(result.height());
        if (basis != nullptr) {
            (*basis) *= Q;
        }
    }
    return result;
}

namespace details {
inline long double get_Wilkinson_shift(const Matrix<long double>& A) {
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

Matrix<long double> apply_qr_for_bidiagonal(const Matrix<long double>&, Matrix<long double>*, Matrix<long double>*,
                                            const long double);

inline long double split(Matrix<long double>& A, Matrix<long double>* left_basis, Matrix<long double>* right_basis,
                         const long double eps = constants::DEFAULT_EPSILON) {
    using Matrix = Matrix<long double>;

    for (size_t ind = 0; ind + 1 < A.width(); ++ind) {
        if (A(ind, ind + 1) <= eps) {
            auto [result1, result2] = split_matrix(A, ind, ind);
            Matrix left_basis1 = Matrix::identity(ind + 1);
            Matrix right_basis1 = Matrix::identity(ind + 1);
            Matrix left_basis2 = Matrix::identity(A.height() - ind - 1);
            Matrix right_basis2 = Matrix::identity(A.height() - ind - 1);
            result1 = apply_qr_for_bidiagonal(result1, &left_basis1, &right_basis1, eps);
            result2 = apply_qr_for_bidiagonal(result2, &left_basis2, &right_basis2, eps);
            A = join_matrix(result1, result2);
            if (left_basis != nullptr) {
                (*left_basis) *= join_matrix(left_basis1, left_basis2);
            }
            if (right_basis != nullptr) {
                (*right_basis) = join_matrix(right_basis1, right_basis2) * (*right_basis);
            }
            return true;
        }
    }
    return false;
}

Matrix<long double> apply_qr_for_bidiagonal(const Matrix<long double>& A, Matrix<long double>* left_basis,
                                            Matrix<long double>* right_basis,
                                            const long double eps = constants::DEFAULT_EPSILON) {
    using Matrix = Matrix<long double>;

    assert(is_bidiagonal(A, eps));

    if (A.height() == 1) {
        if (left_basis != nullptr) {
            (*left_basis) = {{1.0}};
        }
        if (right_basis != nullptr) {
            (*right_basis) = {{1.0}};
        }
        if (A(0, 0) <= eps) {
            return {{0.0}};
        }
        return A;
    }

    Matrix result = A;
    while (!is_diagonal(result, eps)) {
        /*if (split(result, left_basis, right_basis, eps)) {
            return result;
        }*/
        // std::cout << "!!!!!\n";
        // std::cout << result << "\n=================================\n";
        Matrix M = transpose(result) * result;
        long double shift = get_Wilkinson_shift(M);
        auto M_with_shift = M - shift * Matrix::identity(M.height());
        auto [Q, R] = get_QR_decomposition(M_with_shift, eps);

        if (right_basis != nullptr) {
            (*right_basis) *= Q;
        }
        auto l_basis = bidiagonalize_with_right_basis(result, Q, eps);
        if (left_basis != nullptr) {
            (*left_basis) *= l_basis;
        }
        result = transpose(l_basis) * result * Q;
    }
    set_low_values_zero(result);
    return result;
}
}  // namespace details
}  // namespace svd_computation
