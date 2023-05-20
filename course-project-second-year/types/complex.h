#pragma once

#include <iostream>

namespace svd_computation {
class Complex {
   public:
    using Type = long double;

    Complex() = default;
    Complex(int);
    Complex(long long);
    Complex(double);
    Complex(long double);

    Complex(int, int);
    Complex(long long, long long);
    Complex(double, double);
    Complex(long double, long double);

    static Complex exp_form(long double, long double) noexcept;

    Type Re() const noexcept;
    Type Im() const noexcept;

    friend Complex operator+(const Complex&, const Complex&) noexcept;
    Complex& operator+=(const Complex&) noexcept;
    friend Complex operator-(const Complex&, const Complex&) noexcept;
    Complex& operator-=(const Complex&) noexcept;
    friend Complex operator*(const Complex&, const Complex&) noexcept;
    Complex& operator*=(const Complex&) noexcept;
    friend Complex operator/(const Complex&, const Complex&) noexcept;
    Complex& operator/=(const Complex&) noexcept;
    Complex operator-() const noexcept;
    friend bool operator==(const Complex&, const Complex&) noexcept;
    friend bool operator!=(const Complex&, const Complex&) noexcept;

    Type abs() const noexcept;
    Type arg() const noexcept;
    Complex& conjugate() noexcept;

    friend std::istream& operator>>(std::istream&, Complex&) noexcept;
    friend std::ostream& operator<<(std::ostream&, const Complex&) noexcept;

   private:
    Type Re_ = 0.;
    Type Im_ = 0.;
};

Complex::Type abs(const Complex&);
Complex::Type arg(const Complex&);
Complex conjugate(const Complex&);
Complex sqrt(const Complex&);
long double sqrt(const long double&);

}  // namespace svd_computation