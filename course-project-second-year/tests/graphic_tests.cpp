#include "test_performance.h"
#include "test_precision.h"

void performance_tests() {
    std::vector<size_t> ns = {1, 5, 10, 15, 20, 25, 30, 50, 100, 150, 200, 250, 300, 350, 400};

    std::vector<std::pair<int, int>> time_ms = test_performance(ns);
    std::cout << "my time ms:\n";
    for (size_t i = 0; i < ns.size(); i++) {
        std::cout << time_ms[i].first << '\n';
    }
    std::cout << "Eigen's time ms:\n";
    for (size_t i = 0; i < ns.size(); i++) {
        std::cout << time_ms[i].second << '\n';
    }

    std::cout << "precision compare to eigen\n";
    std::vector<long double> precision_eigen = test_precision_eigen(ns);
    for (size_t i = 0; i < ns.size(); i++) {
        std::cout << precision_eigen[i] << '\n';
    }

    std::cout << "precision compare to given matrix\n";
    std::vector<long double> precision = test_precision(ns);
    for (size_t i = 0; i < ns.size(); i++) {
        std::cout << precision[i] << '\n';
    }
}

void operations_test() {
    std::vector<size_t> ops = {1, 5, 10, 15, 20, 30, 50, 75, 100, 250, 500, 750, 1000};

    std::vector<long double> precision = test_precision_operations(ops);
    std::cout << "precision with different operations number\n";
    for (size_t i = 0; i < ops.size(); i++) {
        std::cout << precision[i] << '\n';
    }
}

void epsilon_test() {
    std::vector<long double> eps = {1,     1e-1,  1e-2,  1e-3,  1e-4,  1e-5,  1e-6,  1e-7,  1e-8,  1e-9,
                                    1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18, 1e-19};

    std::vector<long double> precision = test_precision_epsilons(eps);
    std::cout << "precision with different epsilons\n";
    for (size_t i = 0; i < eps.size(); i++) {
        std::cout << precision[i] << '\n';
    }
}
