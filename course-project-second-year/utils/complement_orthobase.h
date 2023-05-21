#pragma once

#include <vector>

#include "../algorithms/constants.h"
#include "../types/vector.h"

namespace svd_computation {
namespace details {
template <typename Type>
void complement_orthobase(std::vector<Vector<Type>>& A,
                          const long double eps = svd_computation::constants::DEFAULT_EPSILON) {
    if (A.empty()) {
        return;
    }

    A.reserve(A.size() + A[0].size());

    for (size_t ind = 0; ind < A[0].size(); ++ind) {
        A.push_back(Vector<Type>::standart_basis(ind, A[0].size(), A[0].orientation()));
    }

    orthonormalize(A, eps);
}
}  // namespace details
}  // namespace svd_computation
