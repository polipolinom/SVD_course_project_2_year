#include "complex.h"

#include <math.h>

#include <cstdlib>
#include <iomanip>
#include <string>

namespace svd_computation {
Complex::Complex() : _Re(0), _Im(0) {}
Complex::Complex(int x) : _Re(x), _Im(0) {}
Complex::Complex(long long x) : _Re(x), _Im(0) {}
Complex::Complex(double x) : _Re(x), _Im(0) {}
Complex::Complex(long double x) : _Re(x), _Im(0) {}

Complex::Complex(int Re, int Im) : _Re(Re), _Im(Im) {}
Complex::Complex(long long Re, long long Im) : _Re(Re), _Im(Im) {}
Complex::Complex(double Re, double Im) : _Re(Re), _Im(Im) {}
Complex::Complex(long double Re, long double Im) : _Re(Re), _Im(Im) {}

Complex Complex::get_complex_by_exp_form(long double abs, long double arg) {
    return Complex(abs * cos(arg), abs * sin(arg));
}

Complex::Type Complex::Re() const {
    return _Re;
}
Complex::Type Complex::Im() const {
    return _Im;
}

Complex Complex::operator+(const Complex& rhs) const {
    return Complex(_Re + rhs._Re, _Im + rhs._Im);
}

Complex& Complex::operator+=(const Complex& rhs) {
    *this = *this + rhs;
    return *this;
}

Complex Complex::operator-(const Complex& rhs) const {
    return Complex(_Re - rhs._Re, _Im - rhs._Im);
}

Complex& Complex::operator-=(const Complex& rhs) {
    *this = *this - rhs;
    return *this;
}

Complex Complex::operator*(const Complex& rhs) const {
    return Complex(_Re * rhs._Re - _Im * rhs._Im, _Re * rhs._Im + _Im * rhs._Re);
}

Complex& Complex::operator*=(const Complex& rhs) {
    *this = *this * rhs;
    return *this;
}

Complex Complex::operator/(const Complex& rhs) const {
    Type abs_rhs = rhs._Re * rhs._Re + rhs._Im * rhs._Im;

    if (abs_rhs == 0.0l) {
        throw std::invalid_argument("Division by zero");
    }

    return Complex((_Re * rhs._Re + _Im * rhs._Im) / abs_rhs, (_Re * rhs._Im - _Im * rhs._Re) / abs_rhs);
}

Complex& Complex::operator/=(const Complex& rhs) {
    *this = *this / rhs;
    return *this;
}

Complex Complex::operator-() const {
    return Complex(-_Re, -_Im);
}

bool Complex::operator==(const Complex& other) const {
    return (_Re == other._Re && _Im == other._Im);
}

long double abs(const Complex& x) {
    return sqrtf(x.Re() * x.Re() + x.Im() * x.Im());
}

long double arg(const Complex& x) {
    return atan2(x.Im(), x.Re());
}

Complex conjugate(const Complex& x) {
    return Complex(x.Re(), -x.Im());
}

Complex sqrt(const Complex& x) {
    long double abs_x = abs(x);
    long double arg_x = arg(x);
    if (arg_x < 0) {
        arg_x += 2 * M_PI;
    }
    return Complex::get_complex_by_exp_form(sqrtf(abs_x), arg_x / 2.0l);
}

std::istream& operator>>(std::istream& in, Complex& x) {
    std::string buf;
    in >> buf;

    std::string::size_type sz;

    long double num1, num2;

    try {
        num1 = std::stold(buf, &sz);
    } catch (std::invalid_argument& e) {
        throw std::invalid_argument("Incorrect format of complex number " + buf);
    } catch (std::out_of_range& e) {
        throw std::out_of_range("Out of range in complex number " + buf);
    }
    if (sz == buf.size()) {
        x = num1;
        return in;
    }

    if (buf[sz] == '*' || buf[sz] == 'i') {
        x = Complex(0.0l, num1);
        return in;
    }

    try {
        num2 = std::stold(buf.substr(sz), &sz);
    } catch (std::invalid_argument& e) {
        throw std::invalid_argument("Incorrect format of complex number " + buf);
    } catch (std::out_of_range& e) {
        throw std::out_of_range("Out of range in complex number " + buf);
    }

    x = Complex(num1, num2);

    return in;
}

std::ostream& operator<<(std::ostream& out, const Complex& x) {
    out << std::fixed << std::setprecision(6);
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