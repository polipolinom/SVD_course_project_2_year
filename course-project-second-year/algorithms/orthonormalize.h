#pragma once

#include <algorithm>
#include <vector>

#include "../types/vector.h"
#include "constants.h"

namespace svd_computation {

template <typename Type>
void orthonormalize(std::vector<Vector<Type>>& system, const long double eps = constants::DEFAULT_EPSILON) {
    for (size_t ind = 0; ind < system.size(); ++ind) {
        if (abs(system[ind]) <= eps) {
            system[ind] = Vector<Type>(Type(0.0), system[ind].size());
            continue;
        }
        system[ind] /= abs(system[ind]);
        for (size_t k = ind + 1; k < system.size(); ++k) {
            assert(system[k].orientation() == system[ind].orientation());
            if (system[k].orientation() == Vector<Type>::Orientation::Vertical) {
                system[k] -= dot_product<Type>(transpose(system[ind]), system[k]) * system[ind];
            } else {
                system[k] -= dot_product<Type>(system[ind], transpose(system[k])) * system[ind];
            }
        }
    }

    system.erase(std::remove(system.begin(), system.end(), Vector<Type>(Type(0.0), system[0].size())), system.end());
}
}  // namespace svd_computation
