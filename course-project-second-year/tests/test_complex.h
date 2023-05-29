#include <gtest/gtest.h>
#include <math.h>

#include "../types/complex.h"
#include "../utils/parser.h"

using namespace svd_computation;

TEST(ComplexTest, ConstructorFromInt) {
    Complex a = Complex(2);
    EXPECT_DOUBLE_EQ(a.Re(), 2.0);
    EXPECT_DOUBLE_EQ(a.Im(), 0.0);
}

TEST(ComplexTest, ConstructorFromLongDouble) {
    Complex a = Complex(2.5);
    EXPECT_DOUBLE_EQ(a.Re(), 2.5);
    EXPECT_DOUBLE_EQ(a.Im(), 0.0);
}

TEST(ComplexTest, ConstructorFromReIm) {
    Complex a = Complex(2.5, -4.5);
    EXPECT_DOUBLE_EQ(a.Re(), 2.5);
    EXPECT_DOUBLE_EQ(a.Im(), -4.5);
}

TEST(ComplexTest, ConstructorExpForm) {
    Complex a = Complex::exp_form(2, 0);
    EXPECT_DOUBLE_EQ(a.Re(), 2.0);
    EXPECT_DOUBLE_EQ(a.Im(), 0);

    Complex b = Complex::exp_form(2, acos(-1.0) / 2.0);
    EXPECT_NEAR(b.Re(), 0.0, 1e-15);
    EXPECT_NEAR(b.Im(), 2.0, 1e-15);
}

TEST(ComplexTest, Parser) {
    Complex a = details::Parser().parse("2+3i");
    EXPECT_DOUBLE_EQ(a.Re(), 2.0);
    EXPECT_DOUBLE_EQ(a.Im(), 3.0);

    Complex a1 = details::Parser().parse("2-3i");
    EXPECT_DOUBLE_EQ(a1.Re(), 2.0);
    EXPECT_DOUBLE_EQ(a1.Im(), -3.0);

    Complex a2 = details::Parser().parse("-2-3i");
    EXPECT_DOUBLE_EQ(a2.Re(), -2.0);
    EXPECT_DOUBLE_EQ(a2.Im(), -3.0);

    Complex b = details::Parser().parse("2");
    EXPECT_DOUBLE_EQ(b.Re(), 2.0);
    EXPECT_DOUBLE_EQ(b.Im(), 0.0);

    Complex b1 = details::Parser().parse("-2");
    EXPECT_DOUBLE_EQ(b1.Re(), -2.0);
    EXPECT_DOUBLE_EQ(b1.Im(), 0.0);

    Complex b2 = details::Parser().parse("+2");
    EXPECT_DOUBLE_EQ(b2.Re(), 2.0);
    EXPECT_DOUBLE_EQ(b2.Im(), 0.0);

    Complex c = details::Parser().parse("3i");
    EXPECT_DOUBLE_EQ(c.Re(), 0.0);
    EXPECT_DOUBLE_EQ(c.Im(), 3.0);

    Complex c1 = details::Parser().parse("+3i");
    EXPECT_DOUBLE_EQ(c1.Re(), 0.0);
    EXPECT_DOUBLE_EQ(c1.Im(), 3.0);

    Complex c2 = details::Parser().parse("-3i");
    EXPECT_DOUBLE_EQ(c2.Re(), 0.0);
    EXPECT_DOUBLE_EQ(c2.Im(), -3.0);
}

TEST(ComplexTest, Sum) {
    Complex a = Complex(2.1, -3.5);
    Complex b = Complex(4.6, 7.8);

    EXPECT_DOUBLE_EQ((a + b).Re(), 6.7);
    EXPECT_DOUBLE_EQ((a + b).Im(), 4.3);

    EXPECT_DOUBLE_EQ((b + a).Re(), 6.7);
    EXPECT_DOUBLE_EQ((b + a).Im(), 4.3);
}

TEST(ComplexTest, Diff) {
    Complex a = Complex(2.1, -3.5);
    Complex b = Complex(4.6, 7.8);

    EXPECT_DOUBLE_EQ((a - b).Re(), -2.5);
    EXPECT_DOUBLE_EQ((a - b).Im(), -11.3);

    EXPECT_DOUBLE_EQ((b - a).Re(), 2.5);
    EXPECT_DOUBLE_EQ((b - a).Im(), 11.3);
}

TEST(ComplexTest, Mult) {
    Complex a = Complex(2.1, 3.5);
    Complex b = Complex(4.6, -1.1);

    EXPECT_DOUBLE_EQ((a * b).Re(), 13.51);
    EXPECT_DOUBLE_EQ((a * b).Im(), 13.79);

    EXPECT_DOUBLE_EQ((b * a).Re(), 13.51);
    EXPECT_DOUBLE_EQ((b * a).Im(), 13.79);
}

TEST(ComplexTest, Div) {
    Complex a = Complex(2.1, 3.5);
    Complex b = Complex(3.0, -4.0);

    EXPECT_DOUBLE_EQ((a / b).Re(), -0.308);
    EXPECT_DOUBLE_EQ((a / b).Im(), 0.756);
}

TEST(ComplexTest, Abs) {
    Complex a = Complex(3.0, -4.0);

    EXPECT_DOUBLE_EQ(abs(a), 5.0);
    EXPECT_DOUBLE_EQ(a.abs(), 5.0);
}

TEST(ComplexTest, Arg) {
    Complex a = Complex(2.0, 0.0);

    EXPECT_NEAR(arg(a), 0.0, 1e-15);
    EXPECT_NEAR(a.arg(), 0.0, 1e-15);

    Complex b = Complex(-2.0, 0.0);

    EXPECT_NEAR(arg(b), acos(-1), 1e-15);
    EXPECT_NEAR(b.arg(), acos(-1), 1e-15);

    Complex c = Complex(0.0, 5.0);

    EXPECT_NEAR(arg(c), acos(-1) / 2.0, 1e-15);
    EXPECT_NEAR(c.arg(), acos(-1) / 2.0, 1e-15);

    Complex d = Complex(0.0, -5.0);

    EXPECT_NEAR(arg(d), -acos(-1) / 2.0, 1e-15);
    EXPECT_NEAR(d.arg(), -acos(-1) / 2.0, 1e-15);
}