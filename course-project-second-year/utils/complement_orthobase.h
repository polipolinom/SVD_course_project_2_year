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

    Vector<Type> zero(A[0].size(), A[0].orientation());

    for (size_t ind = 0; ind < zero.size(); ++ind) {
        zero[ind] = Type(1);
        A.emplace_back(zero);
        zero[ind] = Type(0);
    }

    orthonormalize(A, eps);
}
}  // namespace details
}  // namespace svd_computation