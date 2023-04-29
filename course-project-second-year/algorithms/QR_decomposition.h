#pragma once

#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "orthonormalize.h"

namespace svd_computation {

template <typename Type>
std::pair<Matrix<Type>, Matrix<long double>> get_reduce_QR_decomposition(Matrix<Type>& A,
                                                                         long double eps = DEFAULT_EPSILON) {
    std::vector<Vector<Type>> Q;
    Matrix<long double> R(A.height(), A.width());

    for (size_t ind = 0; ind < A.width(); ++ind) {
        Q.emplace_back(A.column(ind));
        for (size_t k = 0; k < ind; ++k) {
            if (R(k, k) == 0.0) {
                continue;
            }
            R(k, ind) = dot_product<Type>(Q[k].transpose(), Q[ind]);
            Q[ind] -= R(k, ind) * Q[k];
        }
        if (abs<Type>(Q[ind]) <= eps) {
            R(ind, ind) = 0.0;
            Q[ind] = {};
            continue;
        }
        R(ind, ind) = abs<Type>(Q[ind]);
        Q[ind] /= R(ind, ind);
    }

    Q.erase(std::remove_if(Q.begin(), Q.end(), [](Vector<Type>& elem) { return elem.empty(); }), Q.end());

    Matrix<long double> R_result(Q.size(), A.width());

    size_t last_row = 0;

    for (size_t ind = 0; ind < A.height(); ++ind) {
        if (R(ind, ind) == 0.0) {
            continue;
        }
        for (size_t k = 0; k < A.width(); ++k) {
            R_result(last_row, k) = R(ind, k);
        }
        last_row++;
    }

    return {Matrix<Type>::from_vectors(Q), R_result};
}

template <typename Type>
std::pair<Matrix<Type>, Matrix<long double>> get_QR_decomposition(Matrix<Type>& A, long double eps = DEFAULT_EPSILON) {
    auto [Q1, R1] = get_reduce_QR_decomposition(A, eps);

    std::vector<Vector<Type>> Q;
    for (size_t ind = 0; ind < Q1.width(); ++ind) {
        Q.push_back(Q1.column(ind));
    }

    Vector<Type> zero(A.height());

    for (size_t ind = 0; ind < A.height(); ++ind) {
        zero[ind] = Type(1);
        Q.emplace_back(zero);
        zero[ind] = Type(0);
    }

    orthonormalize_vectors(Q, eps);

    Matrix<long double> R(A.height(), A.width());

    for (size_t i = 0; i < R1.height(); ++i) {
        for (size_t j = 0; j < R1.width(); ++j) {
            R(i, j) = R1(i, j);
        }
    }

    return {Matrix<Type>::from_vectors(Q), R};
}

}  // namespace svd_computation