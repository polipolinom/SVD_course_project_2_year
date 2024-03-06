#pragma once

#include <gtest/gtest.h>

#include "../types/complex.h"
#include "../types/vector.h"

using namespace svd_computation;

TEST(VectorTest, Orientation) {
    Vector<int> v = {1, 2, 3};
    EXPECT_EQ(v.orientation(), Vector<int>::Orientation::Vertical);
    EXPECT_NE(transpose(v).orientation(), Vector<int>::Orientation::Vertical);

    transpose(v);
    EXPECT_EQ(v.orientation(), Vector<int>::Orientation::Vertical);

    v.transpose();
    EXPECT_EQ(v.orientation(), Vector<int>::Orientation::Horizontal);
    EXPECT_NE(transpose(v).orientation(), Vector<int>::Orientation::Horizontal);

    v.transpose();
    EXPECT_EQ(v.orientation(), Vector<int>::Orientation::Vertical);
}

TEST(VectorTest, DotProduct) {
    Vector<long double> v1({1, 2, 3, 4}, Vector<long double>::Orientation::Horizontal);
    Vector<long double> v2 = {5, 6, 7, 8};

    EXPECT_DOUBLE_EQ(dot_product(v1, v2), 70.0);
    EXPECT_DOUBLE_EQ(dot_product(transpose(v2), transpose(v1)), 70.0);

    Vector<Complex> u1({Complex(1, 1), Complex(2, 0), Complex(0, 3), Complex(1.2, 5.6)},
                       Vector<Complex>::Orientation::Horizontal);
    Vector<Complex> u2 = {Complex(5, 5), Complex(0.0, -0.1), Complex(3.7, 0.0), Complex(1.2, 5.6)};

    EXPECT_DOUBLE_EQ(dot_product(u1, u2).Re(), 42.8);
    EXPECT_DOUBLE_EQ(dot_product(u1, u2).Im(), -11.3);

    EXPECT_DOUBLE_EQ(dot_product(transpose(u2), transpose(u1)).Re(), 42.8);
    EXPECT_DOUBLE_EQ(dot_product(transpose(u2), transpose(u1)).Im(), 11.3);
}
