#pragma once

#include <vector>

#include "../types/matrix.h"
#include "../utils/complement_orthobase.h"
#include "QR_algorithm.h"
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

void sort_singular_values(Matrix<long double>& sigma, Matrix<long double>& left_basis,
                          Matrix<long double>& right_basis) {
    for (size_t i = 0; i < sigma.height(); ++i) {
        for (size_t j = 0; j < sigma.height() - i - 1; ++j) {
            if (sigma(j, j) < sigma(j + 1, j + 1)) {
                std::swap(sigma(j, j), sigma(j + 1, j + 1));
                swap_columns(left_basis, j, j + 1);
                swap_columns(right_basis, j, j + 1);
            }
        }
    }
}
}  // namespace details

template <typename Type>
std::vector<long double> compute_svd(const Matrix<Type>& A, Matrix<Type>* left_basis = nullptr,
                                     Matrix<Type>* right_basis = nullptr,
                                     const long double eps = constants::DEFAULT_EPSILON) {
    if (A.width() > A.height()) {
        auto sigma = compute_svd(conjugate(A), right_basis, left_basis, eps);

        return sigma;
    }

    Matrix<Type> left_bidiag, right_bidiag;
    auto B = bidiagonalize(A, &left_bidiag, &right_bidiag, eps);

    size_t min_size = B.width();
    Matrix<long double> B1(min_size, min_size);
    for (size_t i = 0; i < min_size; ++i) {
        for (size_t j = 0; j < min_size; ++j) {
            B1(i, j) = B(i, j);
        }
    }

    // std::cout << B1 << "\n============\n";

    Matrix<long double> left_qr = Matrix<long double>::identity(B1.height());
    Matrix<long double> right_qr = Matrix<long double>::identity(B1.width());
    auto result = details::apply_qr_for_bidiagonal(B1, &left_qr, &right_qr, eps);

    for (size_t i = 0; i < B1.width(); ++i) {
        if (result(i, i) < 0.0) {
            for (size_t j = 0; j < B1.height(); ++j) {
                left_qr(j, i) *= -1;
            }
            result(i, i) *= -1;
        }
    }

    details::sort_singular_values(result, left_qr, right_qr);

    // std::cout << result << "\n\n";

    Matrix<long double> new_left_qr(A.height(), A.height());

    for (size_t row = 0; row < A.width(); ++row) {
        for (size_t column = 0; column < A.width(); ++column) {
            new_left_qr(row, column) = left_qr(row, column);
        }
    }

    for (size_t ind = A.height(); ind < A.width(); ++ind) {
        new_left_qr(ind, ind) = 1.0;
    }

    std::vector<long double> sigma;
    for (size_t i = 0; i < result.width(); ++i) {
        sigma.push_back(result(i, i));
    }

    if (left_basis != nullptr) {
        (*left_basis) = left_bidiag * new_left_qr;
    }
    if (right_basis != nullptr) {
        (*right_basis) = right_bidiag * right_qr;
    }

    return sigma;
}
}  // namespace svd_computation
