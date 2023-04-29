#include <iostream>

#include "algorithms/orthonormalize.h"
#include "types/vector.h"

int main() {
    using namespace svd_computation;
    Vector<long double> a = {2, 3};
    Vector<long double> b = {3, 2};
    std::vector<Vector<long double>> v = {a, b};

    orthonormalize_vectors(v);

    for (auto i : v) {
        std::cout << i << std::endl << std::endl;
    }
}