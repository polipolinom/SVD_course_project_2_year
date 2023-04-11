#include <iostream>

namespace svd_computation {
class Complex {
    using Type = long double;

   public:
    Complex();
    Complex(int);
    Complex(long long);
    Complex(double);
    Complex(long double);

    Complex(int, int);
    Complex(long long, long long);
    Complex(double, double);
    Complex(long double, long double);

    static Complex get_complex_by_exp_form(long double, long double);

    Type Re() const;
    Type Im() const;

    Complex operator+(const Complex&) const;
    Complex& operator+=(const Complex&);
    Complex operator-(const Complex&) const;
    Complex& operator-=(const Complex&);
    Complex operator*(const Complex&) const;
    Complex& operator*=(const Complex&);
    Complex operator/(const Complex&) const;
    Complex& operator/=(const Complex&);
    Complex operator-() const;
    bool operator==(const Complex&) const;

   private:
    Type _Re;
    Type _Im;
};

long double abs(const Complex&);
long double arg(const Complex&);
Complex conjugate(const Complex&);
Complex sqrt(const Complex&);

std::istream& operator>>(std::istream&, Complex&);
std::ostream& operator<<(std::ostream&, const Complex&);

}  // namespace svd_computation