#pragma once

#include <utility>
#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/complement_orthobase.h"
#include "constants.h"
#include "householder_reflections.h"
#include "orthonormalize.h"

namespace svd_computation {

template <typename Type>
std::pair<Matrix<Type>, Matrix<Type>> get_QR_decomposition(const Matrix<Type>& A,
                                                           const long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> R = A;
    Matrix<Type> Q = Matrix<Type>::identity(A.height());
    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        if (ind + 1 == A.height()) {
            break;
        }
        left_reflection(R, ind, ind, &Q, eps);
    }
    Q.conjugate();
    return {Q, R};
}

}  // namespace svd_computation
