#pragma once

#include <algorithm>
#include <vector>

#include "../types/vector.h"
#include "constants.h"

namespace svd_computation {

template <typename Type>
void orthonormalize(std::vector<Vector<Type>>& system, const long double eps = constants::DEFAULT_EPSILON) {
    for (size_t ind = 0; ind < system.size(); ++ind) {
        if (abs<Type>(system[ind]) <= eps) {
            system[ind] = {};
            continue;
        }
        system[ind] /= abs<Type>(system[ind]);
        for (size_t k = ind + 1; k < system.size(); ++k) {
            assert(system[k].is_vertical() == system[ind].is_vertical());
            system[k] -= dot_product<Type>(system[ind].transpose(), system[k]) * system[ind];
        }
    }

    system.erase(std::remove_if(system.begin(), system.end(), [](Vector<Type>& elem) { return elem.empty(); }),
                 system.end());
}
}  // namespace svd_computation