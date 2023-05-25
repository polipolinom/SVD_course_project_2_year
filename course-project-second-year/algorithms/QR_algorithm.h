#pragma once

#include <cassert>

#include "../types/matrix.h"
#include "../utils/checking_matrices.h"
#include "../utils/matrix_split_join.h"
#include "../utils/set_values.h"
#include "QR_decomposition.h"
#include "bidiagonalization.h"
#include "constants.h"
#include "givens_rotation.h"

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
inline long double get_Wilkinson_shift(const long double n_1n_1, const long double n_1n, const long double nn) {
    long double shift = nn;
    long double sigma = (n_1n_1 - nn) / 2.0;

    long double coef = -1.0;

    if (sigma < 0) {
        coef = 1.0;
    }
    shift += coef * (n_1n * n_1n) / (abs(sigma) + sqrtl(sigma * sigma + n_1n * n_1n));

    return shift;
}

Matrix<long double> apply_qr_for_bidiagonal(const Matrix<long double>&, Matrix<long double>*, Matrix<long double>*,
                                            const long double);

inline long double split(Matrix<long double>& A, Matrix<long double>* left_basis, Matrix<long double>* right_basis,
                         const long double eps = constants::DEFAULT_EPSILON) {
    // std::cout << "split\n";

    using Matrix = Matrix<long double>;

    for (size_t ind = 0; ind + 1 < A.width(); ++ind) {
        if (abs(A(ind, ind + 1)) < eps) {
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
                (*right_basis) *= join_matrix(right_basis1, right_basis2);
            }
            return true;
        }
    }
    return false;
}

inline void multiplyRightGivens(Matrix<long double>& A, long double c, long double s, int i, int j) {
    for (size_t k = 0; k < A.height(); ++k) {
        long double a = A(k, i);
        long double b = A(k, j);

        A(k, i) = c * a - s * b;
        A(k, j) = s * a + c * b;
    }
}

inline void multiplyLeftGivens(Matrix<long double>& A, long double c, long double s, int i, int j) {
    for (size_t k = 0; k < A.width(); ++k) {
        long double a = A(i, k);
        long double b = A(j, k);

        A(i, k) = c * a - s * b;
        A(j, k) = s * a + c * b;
    }
}

Matrix<long double> apply_qr_for_bidiagonal(const Matrix<long double>& A, Matrix<long double>* left_basis,
                                            Matrix<long double>* right_basis,
                                            const long double eps = constants::DEFAULT_EPSILON) {
    using Matrix = Matrix<long double>;

    assert(is_bidiagonal(A, eps));

    if (A.height() == 1) {
        if (left_basis != nullptr) {
            (*left_basis) = {{1.0}};
            if (A(0, 0) < -eps) {
                (*left_basis) = {{-1.0}};
            }
        }
        if (right_basis != nullptr) {
            (*right_basis) = {{1.0}};
        }
        if (abs(A(0, 0)) <= eps) {
            return {{0.0}};
        }
        if (A(0, 0) < -eps) {
            return {{-A(0, 0)}};
        }
        return A;
    }

    Matrix result = A;
    const int n = result.height();
    int operations = 0;
    while (!is_diagonal(result, eps) && operations < 50 * result.height()) {
        operations++;
        // std::cout << result << "\n!!!!!!!!!!!!!!!\n";
        if (split(result, left_basis, right_basis, eps)) {
            return result;
        }

        for (size_t ind = 0; ind + 1 < result.height(); ++ind) {
            if (abs(result(ind, ind)) <= eps) {
                for (size_t k = ind + 1; k < result.height(); ++k) {
                    auto [cos, sin] = get_givens_rotation(result(k, k), result(ind, k), eps);

                    multiplyLeftGivens(result, cos, -sin, ind, k);
                    if (left_basis != nullptr) {
                        multiplyRightGivens(*left_basis, cos, -sin, ind, k);
                    }
                }
                split(result, left_basis, right_basis, eps);
                return result;
            }
        }

        long double n_2n_1 = 0;
        if (n - 3 >= 0) {
            n_2n_1 = result(n - 3, n - 2);
        }
        long double n_1n_1 = result(n - 2, n - 2);
        long double n_1n = result(n - 2, n - 1);
        long double nn = result(n - 1, n - 1);

        long double shift =
            get_Wilkinson_shift(n_2n_1 * n_2n_1 + n_1n_1 * n_1n_1, n_1n_1 * n_1n, nn * nn + n_1n * n_1n);

        for (size_t ind = 0; ind + 1 < result.height(); ++ind) {
            if (ind == 0) {
                auto [cos, sin] =
                    get_givens_rotation(result(0, 0) * result(0, 0) - shift, result(0, 0) * result(0, 1), eps);

                multiplyRightGivens(result, cos, sin, 0, 1);
                if (right_basis != nullptr) {
                    multiplyRightGivens(*right_basis, cos, sin, 0, 1);
                }
            } else {
                auto [cos, sin] = get_givens_rotation(result(ind - 1, ind), result(ind - 1, ind + 1), eps);

                multiplyRightGivens(result, cos, sin, ind, ind + 1);
                if (right_basis != nullptr) {
                    multiplyRightGivens(*right_basis, cos, sin, ind, ind + 1);
                }
            }
            auto [cos, sin] = get_givens_rotation(result(ind, ind), result(ind + 1, ind));

            multiplyLeftGivens(result, cos, sin, ind, ind + 1);
            if (left_basis != nullptr) {
                multiplyRightGivens(*left_basis, cos, sin, ind, ind + 1);
            }
        }

        // std::cout << result << "\n##################\n";
    }
    set_low_values_zero(result);
    return result;
}
}  // namespace details
}  // namespace svd_computation
