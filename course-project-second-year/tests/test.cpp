#include <gtest/gtest.h>

#include "../types/complex.h"

using namespace svd_computation;

TEST(ComplexTest, ConstructorFromIntWorks) {
    Complex a = Complex(2);
    EXPECT_EQ(a.Re(), 2);
}

TEST(ComplexTest, MalformedTest) {
    Complex a = Complex(2);
    EXPECT_EQ(a.Re(), 3);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}