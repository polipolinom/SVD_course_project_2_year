#pragma once

#include <vector>

#include "vector_and_matrix.h"

namespace svd_computation {
template <typename Type>
class Vector {
   public:
    Vector() {}

    Vector(size_t length) : _data(std::vector<Type>(length, 0)) {}

    Vector(Type lambda, size_t length) : _data(std::vector<Type>(length, lambda)) {}

    Vector(std::vector<Type> data) : _data(data) {}

    size_t size() {
        return _data.size();
    }

    Type& operator[](size_t ind) {
        return _data[ind];
    }

    const Type& operator[](size_t ind) const {
        return (const Type&)_data[ind];
    }

    Vector operator+(const Vector& rhs) const {
        if (rhs.size() != _data.size()) {
            throw std::invalid_argument("Incorrected sizes in vector summation");
        }
        Vector res = *this;
        for (size_t i = 0; i < res.size(); i++) {
            res[i] += rhs[i];
        }
        return res;
    }

    Vector& operator+=(const Vector& rhs) const {
        *this = *this + rhs;
        return *this;
    }

    Vector operator-(const Vector& rhs) const {
        if (rhs.size() != _data.size()) {
            throw std::invalid_argument("Incorrected sizes in vector substraction");
        }
        Vector res = *this;
        for (size_t i = 0; i < res.size(); i++) {
            res[i] -= rhs[i];
        }
        return res;
    }

    Vector& operator-=(const Vector& rhs) const {
        *this = *this - rhs;
        return *this;
    }

    // scalar product
    Type operator*(const Vector& rhs) const {
        if (rhs.size() != _data.size()) {
            throw std::invalid_argument("Incorrected sizes in scalar product");
        }
        Type res = 0;
        for (size_t i = 0; i < res.size(); i++) {
            res += _data[i] * rhs[i];
        }
        return res;
    }

    Vector& operator*=(const Type& rhs) const {
        *this = rhs * (*this);
        return *this;
    }

    Vector& operator/=(const Type& rhs) const {
        *this = (*this) / rhs;
        return *this;
    }

    Matrix<Type> T() {
        std::vector<std::vector<Type>> res(1, std::vector<Type>(_data.size(), 0));
        for (size_t ind = 0; ind < _data.size(); ind++) {
            res[0][ind] = _data[ind];
        }
        return res;
    }

   private:
    std::vector<Type> _data;
};

template <typename Type>
Vector<Type> operator*(const Type& lambda, const Vector<Type>& v) {
    Vector<Type> res = v;
    for (auto& num : res) {
        num *= lambda;
    }
    return res;
}

template <typename Type>
Vector<Type> operator/(const Vector<Type>& v, const Type& lambda) {
    Vector<Type> res = v;
    for (auto& num : res) {
        num /= lambda;
    }
    return res;
}

template <typename Type>
std::istream& operator>>(std::istream& in, Vector<Type>& v) {
    size_t sz;
    in >> sz;
    v = Vector(sz);
    for (size_t ind = 0; ind < sz; ind++) {
        in >> v[ind];
    }
    return in;
}

template <typename Type>
std::ostream& operator<<(std::ostream& out, const Vecto<rType>& v) {
    for (size_t ind = 0; ind < v.size(); ind++) {
        out << v[ind];
        if (ind + 1 < v.size()) {
            out << "\n";
        }
    }
    return out;
}
