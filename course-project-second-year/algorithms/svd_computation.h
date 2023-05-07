#pragma once

#include <vector>

#include "../types/matrix.h"
#include "../utils/complement_orthobase.h"
#include "QR_decomposition.h"
#include "bidiagonalization.h"
#include "constants.h"

namespace svd_computation {
namespace details {
void swap_columns(Matrix<long double>& A, int ind1, int ind2) {
    assert(ind1 >= 0 && ind1 < A.width());
    assert(ind2 >= 0 && ind2 < A.width());
    for (int i = 0; i < A.height(); ++i) {
        std::swap(A(i, ind1), A(i, ind2));
    }
}

void sort_singular_values(Matrix<long double>& sigma, Matrix<long double>& basis) {
    for (size_t i = 0; i < sigma.height(); ++i) {
        for (size_t j = 0; j < sigma.height() - i - 1; ++j) {
            if (sigma(j, j) < sigma(j + 1, j + 1)) {
                std::swap(sigma(j, j), sigma(j + 1, j + 1));
                swap_columns(basis, j, j + 1);
            }
        }
    }
}
}  // namespace details

template <typename Type>
Matrix<long double> compute_svd(const Matrix<Type>& A, Matrix<Type>* left_basis = nullptr,
                                Matrix<Type>* right_basis = nullptr,
                                const long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> left_bidiag, right_bidiag;
    auto B = bidiagonalize(A, &left_bidiag, &right_bidiag, eps);

    size_t min_size = std::min(B.height(), B.width());
    Matrix<long double> B1(min_size, min_size);
    for (size_t i = 0; i < min_size; ++i) {
        B1(i, i) = B(i, i);
        if (i + 1 < min_size) {
            B1(i, i + 1) = B(i, i + 1);
        }
    }

    if (A.width() < A.height()) {
        /*        auto M = B1.transpose() * B1;
                Matrix<long double> right_qr = Matrix<long double>::ones(M.width());
                auto result = details::apply_qr_for_tridiagonal(M, &right_qr, eps);
                details::sort_singular_values(result, right_qr);
                if (right_basis != nullptr) {
                    (*right_basis) = right_qr.transpose() * right_bidiag;
                }
                Matrix<long double> sigma(A.height(), A.width());
                for (size_t i = 0; i < result.height(); ++i) {
                    sigma(i, i) = sqrtl(result(i, i));
                }
                std::vector<Vector<long double>> left_qr;
                for (size_t i = 0; i < result_height(); ++i) {
                    if (sigma(i, i) == 0.0) {
                        break;
                    }
                    left_qr.emplace_back((B1 * right_qr.transpose()).column(i) / sigma(i, i));
                }
                complement_orthobase(left_qr);
                if (left_basis != nullptr) {
                    (*left_basis) = left_bidiag * Matrix<long double>::from_vectors(left_qr);
                }
                return sigma;
        */
        B = B.transpose();
        B1 = B1.transpose();
    }

    auto M = B1 * B1.transpose();
    Matrix<long double> left_qr = Matrix<long double>::ones(M.height());
    auto result = details::apply_qr_for_tridiagonal(M, &left_qr, eps);

    details::sort_singular_values(result, left_qr);

    Matrix<long double> sigma(B.height(), B.width());
    for (size_t i = 0; i < result.height(); ++i) {
        sigma(i, i) = sqrtl(result(i, i));
    }

    std::vector<Vector<long double>> right_qr;
    for (size_t i = 0; i < result.height(); ++i) {
        if (sigma(i, i) == 0.0) {
            break;
        }
        right_qr.emplace_back((B.transpose() * left_qr).column(i) / sigma(i, i));
    }
    details::complement_orthobase(right_qr, eps);

    if (A.width() < A.height()) {
        sigma = sigma.transpose();

        if (left_basis != nullptr) {
            (*left_basis) = left_bidiag * Matrix<long double>::from_vectors(right_qr).transpose();
        }

        if (right_basis != nullptr) {
            (*right_basis) = right_bidiag * left_qr.transpose();
        }
    } else {
        if (left_basis != nullptr) {
            (*left_basis) = left_bidiag * left_qr;
        }
        if (right_basis != nullptr) {
            (*right_basis) = right_bidiag * Matrix<long double>::from_vectors(right_qr);
        }
    }

    return sigma;
}
}  // namespace svd_computation