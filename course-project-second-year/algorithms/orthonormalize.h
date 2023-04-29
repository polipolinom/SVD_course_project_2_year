#pragma once

#include <algorithm>
#include <vector>

#include "../types/vector.h"

namespace svd_computation {

constexpr long double DEFAULT_EPSILON = 1e-6;

template <typename Type>
void orthonormalize_vectors(std::vector<Vector<Type>>& system, long double eps = DEFAULT_EPSILON) {
    for (size_t ind = 0; ind < system.size(); ++ind) {
        for (size_t k = 0; k < ind; ++k) {
            if (system[k].empty()) {
                continue;
            }
            system[ind] -= dot_product<Type>(system[k].transpose(), system[ind]) * system[k];
        }
        if (abs<Type>(system[ind]) <= eps) {
            system[ind] = {};
            continue;
        }
        system[ind] /= abs<Type>(system[ind]);
    }

    system.erase(std::remove_if(system.begin(), system.end(), [](Vector<Type>& elem) { return elem.empty(); }),
                 system.end());
}
}  // namespace svd_computation