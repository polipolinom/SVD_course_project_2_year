#include <iostream>
#include <vector>

#include "../types/vector.h"

namespace svd_computation {

constexpr long double DEFAULT_EPSILON = 1e-6;

template <typename Type>
void orthonormalize_vectors(std::vector<Vector<Type>>& system, long double eps = DEFAULT_EPSILON) {
    if (system.empty()) {
        return;
    }

    for (size_t ind = 0; ind < system.size(); ++ind) {
        for (size_t k = 0; k < ind; ++k) {
            system[ind] -= dot_product<Type>(system[k].transpose(), system[ind]) * system[k];
        }
        assert(abs(system[ind]) > eps);
        system[ind] /= abs<Type>(system[ind]);
    }
}
}  // namespace svd_computation