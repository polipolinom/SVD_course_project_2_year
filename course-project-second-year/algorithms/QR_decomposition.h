#pragma once

#include <iostream>
#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/complement_orthobase.h"
#include "constants.h"
#include "orthonormalize.h"

namespace svd_computation {

template <typename Type>
std::pair<Matrix<Type>, Matrix<Type>> get_reduce_QR_decomposition(const Matrix<Type>& A,
                                                                  const long double eps = constants::DEFAULT_EPSILON) {
    std::vector<Vector<Type>> Q;
    Matrix<Type> R(A.height(), A.width());

    for (size_t ind = 0; ind < A.width(); ++ind) {
        Q.emplace_back(A.column(ind));
    }

    for (size_t ind = 0; ind < Q.size(); ++ind) {
        if (abs<Type>(Q[ind]) <= eps) {
            R(ind, ind) = Type(0.0);
            Q[ind] = Vector<Type>(Type(0.0), Q[ind].size());
            continue;
        }

        R(ind, ind) = abs(Q[ind]);
        Q[ind] /= R(ind, ind);

        for (size_t k = ind + 1; k < Q.size(); ++k) {
            R(ind, k) = Type(dot_product<Type>(transpose(Q[ind]), Q[k]));
            Q[k] -= R(ind, k) * Q[ind];
        }
    }

    Q.erase(std::remove(Q.begin(), Q.end(), Vector<Type>(Type(0.0), Q[0].size())), Q.end());

    Matrix<Type> R_result(Q.size(), A.width());

    size_t last_row = 0;

    for (size_t ind = 0; ind < A.height(); ++ind) {
        if (R(ind, ind) == Type(0.0)) {
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
std::pair<Matrix<Type>, Matrix<Type>> get_QR_decomposition(const Matrix<Type>& A,
                                                           const long double eps = constants::DEFAULT_EPSILON) {
    auto [Q1, R1] = get_reduce_QR_decomposition(A, eps);

    std::vector<Vector<Type>> Q;
    for (size_t ind = 0; ind < Q1.width(); ++ind) {
        Q.push_back(Q1.column(ind));
    }

    details::complement_orthobase<Type>(Q, eps);

    Matrix<Type> R(A.height(), A.width());

    for (size_t i = 0; i < R1.height(); ++i) {
        for (size_t j = 0; j < R1.width(); ++j) {
            R(i, j) = R1(i, j);
        }
    }

    return {Matrix<Type>::from_vectors(Q), R};
}

}  // namespace svd_computation