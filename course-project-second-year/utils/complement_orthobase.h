#pragma once

#include <random>
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

template <typename Type>
Vector<Type> add_one_vector(const std::vector<Vector<Type>>& A, const size_t length,
                            const typename Vector<Type>::Orientation orientation,
                            const long double eps = svd_computation::constants::DEFAULT_EPSILON) {
    assert(length > A.size());

    std::default_random_engine gen;
    std::uniform_real_distribution<long double> distribution(0.0, 1.0);

    Vector<Type> result(length, orientation);
    for (size_t ind = 0; ind < length; ++ind) {
        /*for (size_t ind = 0; ind < result.size(); ++ind) {
            result[ind] = distribution(gen);
        }*/
        result = Vector<Type>::standart_basis(ind, length, orientation);
        for (auto v : A) {
            if (v.orientation() == Vector<Type>::Orientation::Vertical) {
                result -= dot_product(transpose(v), result) * v;
            } else {
                result -= dot_product(v, transpose(result)) * v;
            }
        }
        if (abs(result) <= eps) {
            continue;
        }
        break;
    }

    return result / abs(result);
}
}  // namespace details
}  // namespace svd_computation
