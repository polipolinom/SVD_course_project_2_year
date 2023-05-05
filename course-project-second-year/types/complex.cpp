#define _USE_MATH_DEFINES

#include "complex.h"

#include <math.h>

#include <cstdlib>
#include <iomanip>
#include <string>

#include "../utils/parser.h"

namespace svd_computation {
Complex::Complex(int x) : Re_(x), Im_(0) {}
Complex::Complex(long long x) : Re_(x), Im_(0) {}
Complex::Complex(double x) : Re_(x), Im_(0) {}
Complex::Complex(long double x) : Re_(x), Im_(0) {}

Complex::Complex(int Re, int Im) : Re_(Re), Im_(Im) {}
Complex::Complex(long long Re, long long Im) : Re_(Re), Im_(Im) {}
Complex::Complex(double Re, double Im) : Re_(Re), Im_(Im) {}
Complex::Complex(long double Re, long double Im) : Re_(Re), Im_(Im) {}

Complex Complex::exp_form(long double abs, long double arg) noexcept {
    return Complex(abs * cos(arg), abs * sin(arg));
}

Complex::Type Complex::Re() const noexcept {
    return Re_;
}
Complex::Type Complex::Im() const noexcept {
    return Im_;
}

Complex& Complex::operator+=(const Complex& rhs) noexcept {
    *this = Complex(Re_ + rhs.Re_, Im_ + rhs.Im_);
    return *this;
}

Complex operator+(const Complex& lhs, const Complex& rhs) noexcept {
    Complex result = lhs;
    result += rhs;
    return result;
}

Complex& Complex::operator-=(const Complex& rhs) noexcept {
    *this = Complex(Re_ - rhs.Re_, Im_ - rhs.Im_);
    return *this;
}

Complex operator-(const Complex& lhs, const Complex& rhs) noexcept {
    Complex result = lhs;
    result -= rhs;
    return result;
}

Complex& Complex::operator*=(const Complex& rhs) noexcept {
    *this = Complex(Re_ * rhs.Re_ - Im_ * rhs.Im_, Re_ * rhs.Im_ + Im_ * rhs.Re_);
    return *this;
}

Complex operator*(const Complex& lhs, const Complex& rhs) noexcept {
    Complex result = lhs;
    result *= rhs;
    return result;
}

Complex& Complex::operator/=(const Complex& rhs) noexcept {
    Type abs_rhs = rhs.Re_ * rhs.Re_ + rhs.Im_ * rhs.Im_;
    *this = Complex((Re_ * rhs.Re_ + Im_ * rhs.Im_) / abs_rhs, (Im_ * rhs.Re_ - Re_ * rhs.Im_) / abs_rhs);
    if (isnan(Re_) || isnan(Im_)) {
        *this = Complex(NAN, NAN);
    }
    return *this;
}

Complex operator/(const Complex& lhs, const Complex& rhs) noexcept {
    Complex result = lhs;
    result /= rhs;
    return result;
}

Complex Complex::operator-() const noexcept {
    return Complex(-Re_, -Im_);
}

bool operator==(const Complex& lhs, const Complex& rhs) noexcept {
    return (lhs.Re_ == rhs.Re_ && lhs.Im_ == rhs.Im_);
}

bool operator!=(const Complex& lhs, const Complex& rhs) noexcept {
    return !(rhs == lhs);
}

long double Complex::abs() const noexcept {
    return sqrtf(Re_ * Re_ + Im_ * Im_);
}

long double abs(const Complex& num) {
    return num.abs();
}

long double Complex::arg() const noexcept {
    return atan2(Im_, Re_);
}

long double arg(const Complex& num) {
    return num.arg();
}

Complex Complex::conjugate() const noexcept {
    return Complex(Re_, -Im_);
}

Complex conjugate(const Complex& num) {
    return num.conjugate();
}

long double conjugate(const long double& num) {
    return num;
}

Complex sqrt(const Complex& num) {
    long double abs_num = abs(num);
    long double arg_num = arg(num);
    if (arg_num < 0) {
        arg_num += 2 * M_PI;
    }
    return Complex::exp_form(sqrtl(abs_num), arg_num / 2.0L);
}

long double sqrt(const long double& num) {
    return sqrtl(num);
}

std::istream& operator>>(std::istream& in, Complex& num) noexcept {
    std::string buf;
    in >> buf;
    num = details::Parser().parse(buf);
    return in;
}

std::ostream& operator<<(std::ostream& out, const Complex& x) noexcept {
    long double eps = 1e-6;
    if (std::fabs(x.Re()) < eps && std::fabs(x.Im()) < eps) {
        out << 0.;
        return out;
    }

    if (std::fabs(x.Re()) >= eps) {
        out << x.Re();
    }
    if (std::fabs(x.Im()) < eps) {
        return out;
    }
    if (x.Im() > 0 && std::fabs(x.Re()) >= eps) {
        out << "+";
    }
    out << x.Im() << "i";
    return out;
}

}  // namespace svd_computation