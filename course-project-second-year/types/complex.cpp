#define _USE_MATH_DEFINES

#include "complex.h"

#include <math.h>

#include <cstdlib>
#include <iomanip>
#include <string>

#include "../utils/parser.h"

namespace svd_computation {
Complex::Complex(int x) : Re_(x) {}
Complex::Complex(long long x) : Re_(x) {}
Complex::Complex(double x) : Re_(x) {}
Complex::Complex(long double x) : Re_(x) {}

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
    Re_ += rhs.Re_;
    Im_ += rhs.Im_;
    return *this;
}

Complex operator+(const Complex& lhs, const Complex& rhs) noexcept {
    Complex result = lhs;
    result += rhs;
    return result;
}

Complex& Complex::operator-=(const Complex& rhs) noexcept {
    Re_ -= rhs.Re_;
    Im_ -= rhs.Im_;
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
    if (isnan(Re_) || isnan(Im_)) {  // I think, +-INF makes more sence than nan
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
    if (isnan(lhs.Re_) || isnan(lhs.Im_)) {
        if (isnan(rhs.Re_) || isnan(rhs.Im_)) {
            return true;
        }
        return false;
    }
    if (isnan(rhs.Re_) || isnan(rhs.Im_)) {
        return false;
    }
    return (lhs.Re_ == rhs.Re_ && lhs.Im_ == rhs.Im_);
}

bool operator!=(const Complex& lhs, const Complex& rhs) noexcept {
    return !(lhs == rhs);
}
Complex::Type Complex::abs() const noexcept {
    return sqrtl(Re_ * Re_ + Im_ * Im_);
}

Complex::Type abs(const Complex& num) {
    return num.abs();
}

Complex::Type Complex::arg() const noexcept {
    return atan2(Im_, Re_);
}

Complex::Type arg(const Complex& num) {
    return num.arg();
}

Complex& Complex::conjugate() noexcept {
    Im_ = -Im_;
    return *this;
}

Complex conjugate(const Complex& num) {
    return {num.Re(), -num.Im()};
}

long double conjugate(const long double& num) {
    return num;
}

// returns sqrt with argument in [0, pi)
Complex sqrt(const Complex& num) {
    Complex::Type abs_num = abs(num);
    Complex::Type arg_num = arg(num);
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
    if (x.Re() == 0.0 && x.Im() == 0.0) {
        out << 0;
        return out;
    }
    if (x.Re() != 0.0) {
        out << x.Re();
    }
    if (x.Im() == 0.0) {
        return out;
    }
    if (x.Im() > 0.0) {
        out << "+";
    }
    out << x.Im() << "i";
    return out;
}

}  // namespace svd_computation
