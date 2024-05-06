#pragma once

#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "../algorithms/orthonormalize.h"
#include "../types/complex.h"
#include "../types/vector.h"
#include "random_objects.h"

using namespace svd_computation;

TEST(OrtonormalizationTest, TestForLongDouble) {
    using Vector = Vector<long double>;

    int operations = 1000;
    int max_sz = 1000;
    int max_cnt = 1000;
    long double max_number = 1e5;

    long double eps = 1e-10;

    std::mt19937 rnd;

    while (operations-- > 0) {
        std::vector<Vector> system(rnd() % max_cnt + 1);
        int sz = rnd() % max_sz + 1;

        for (int i = 0; i < system.size(); ++i) {
            system[i] = get_random_double_vector(sz, sz, max_number);
        }

        orthonormalize(system, eps);

        for (size_t i = 0; i < system.size(); ++i) {
            for (size_t j = 0; j < system.size(); ++j) {
                if (i != j) {
                    EXPECT_NEAR(dot_product(transpose(system[i]), system[j]), 0.0, eps);
                    continue;
                }

                EXPECT_NEAR(abs(system[i]), 1.0, eps);
            }
        }
    }
}

TEST(OrtonormalizationTest, TestForLongDoubleMaxSize) {
    using Vector = Vector<long double>;

    int operations = 300;
    int max_sz = 1000;
    int max_cnt = 1100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    std::mt19937 rnd;

    while (operations-- > 0) {
        std::vector<Vector> system(max_cnt);

        for (int i = 0; i < system.size(); ++i) {
            system[i] = get_random_double_vector(max_sz, max_sz, max_number);
        }

        orthonormalize(system, eps);

        for (size_t i = 0; i < system.size(); ++i) {
            for (size_t j = 0; j < system.size(); ++j) {
                if (i != j) {
                    EXPECT_NEAR(dot_product(transpose(system[i]), system[j]), 0.0, eps);
                    continue;
                }

                EXPECT_NEAR(abs(system[i]), 1.0, eps);
            }
        }
    }
}

TEST(OrtonormalizationTest, TestForLongDoubleLowValues) {
    using Vector = Vector<long double>;

    int operations = 1000;
    int max_sz = 1000;
    int max_cnt = 1000;
    long double max_number = 1e5;

    long double eps = 1e-10;

    std::mt19937 rnd;

    while (operations-- > 0) {
        std::vector<Vector> system(rnd() % max_cnt + 1);
        int sz = rnd() % max_sz + 1;

        for (int i = 0; i < system.size(); ++i) {
            system[i] = get_random_double_vector(sz, sz, 1.0);
        }

        orthonormalize(system, eps);

        for (size_t i = 0; i < system.size(); ++i) {
            for (size_t j = 0; j < system.size(); ++j) {
                if (i != j) {
                    EXPECT_NEAR(dot_product(transpose(system[i]), system[j]), 0.0, eps);
                    continue;
                }

                EXPECT_NEAR(abs(system[i]), 1.0, eps);
            }
        }
    }
}

TEST(OrtonormalizationTest, TestForComplex) {
    using Vector = Vector<Complex>;

    int operations = 500;
    int max_sz = 1000;
    int max_cnt = 1000;
    long double max_number = 1e5;

    long double eps = 1e-10;

    std::mt19937 rnd;

    while (operations-- > 0) {
        std::vector<Vector> system(rnd() % max_cnt + 1);
        int sz = rnd() % max_sz + 1;

        for (int i = 0; i < system.size(); ++i) {
            system[i] = get_random_complex_vector(sz, sz, max_number);
        }

        orthonormalize(system, eps);

        for (size_t i = 0; i < system.size(); ++i) {
            for (size_t j = 0; j < system.size(); ++j) {
                if (i != j) {
                    EXPECT_NEAR(abs(dot_product(transpose(system[i]), system[j])), 0.0, eps);
                    continue;
                }

                EXPECT_NEAR(abs(system[i]), 1.0, eps);
            }
        }
    }
}

TEST(OrtonormalizationTest, TestForComplexMaxSize) {
    using Vector = Vector<Complex>;

    int operations = 150;
    int max_sz = 1000;
    int max_cnt = 1100;
    long double max_number = 1e5;

    long double eps = 1e-10;

    std::mt19937 rnd;

    while (operations-- > 0) {
        std::vector<Vector> system(max_cnt);

        for (int i = 0; i < system.size(); ++i) {
            system[i] = get_random_complex_vector(max_sz, max_sz, max_number);
        }

        orthonormalize(system, eps);

        for (size_t i = 0; i < system.size(); ++i) {
            for (size_t j = 0; j < system.size(); ++j) {
                if (i != j) {
                    EXPECT_NEAR(abs(dot_product(transpose(system[i]), system[j])), 0.0, eps);
                    continue;
                }

                EXPECT_NEAR(abs(system[i]), 1.0, eps);
            }
        }
    }
}

TEST(OrtonormalizationTest, TestForComplexLowValues) {
    using Vector = Vector<Complex>;

    int operations = 500;
    int max_sz = 1000;
    int max_cnt = 1000;
    long double max_number = 1e5;

    long double eps = 1e-10;

    std::mt19937 rnd;

    while (operations-- > 0) {
        std::vector<Vector> system(rnd() % max_cnt + 1);
        int sz = rnd() % max_sz + 1;

        for (int i = 0; i < system.size(); ++i) {
            system[i] = get_random_complex_vector(sz, sz, 1.0);
        }

        orthonormalize(system, eps);

        for (size_t i = 0; i < system.size(); ++i) {
            for (size_t j = 0; j < system.size(); ++j) {
                if (i != j) {
                    EXPECT_NEAR(abs(dot_product(transpose(system[i]), system[j])), 0.0, eps);
                    continue;
                }

                EXPECT_NEAR(abs(system[i]), 1.0, eps);
            }
        }
    }
}
