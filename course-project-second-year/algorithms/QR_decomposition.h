#pragma once

#include <iostream>
#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/complement_orthobase.h"
#include "bidiagonalization.h"
#include "constants.h"
#include "orthonormalize.h"

namespace svd_computation {

template <typename Type>
std::pair<Matrix<Type>, Matrix<Type>> get_QR_decomposition(const Matrix<Type>& A,
                                                           const long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> R = A;
    Matrix<Type> Q = Matrix<Type>::identity(A.height());
    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        details::left_reflection(R, ind, &Q, eps);
    }
    Q.conjugate();
    return {Q, R};
}

}  // namespace svd_computation