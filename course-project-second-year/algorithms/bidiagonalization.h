#pragma once

#include <cmath>
#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/complement_orthobase.h"
#include "../utils/set_values.h"
#include "householder_reflections.h"

namespace svd_computation {
template <typename Type>
Matrix<Type> bidiagonalize(const Matrix<Type>& A, Matrix<Type>* left_basis = nullptr,
                           Matrix<Type>* right_basis = nullptr, const long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> B = A;
    if (left_basis != nullptr) {
        *left_basis = Matrix<Type>::identity(A.height());
    }
    if (right_basis != nullptr) {
        *right_basis = Matrix<Type>::identity(A.width());
    }

    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        left_reflection(B, ind, ind, left_basis, eps);
        if (ind + 1 < A.width()) {
            right_reflection(B, ind, ind + 1, right_basis, eps);
        }
    }

    if (left_basis != nullptr) {
        (*left_basis).conjugate();
    }

    return B;
}
}  // namespace svd_computation
