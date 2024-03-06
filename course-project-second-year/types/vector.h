#pragma once

#include <math.h>

#include <cassert>
#include <vector>

#include "complex.h"

namespace svd_computation {
template <typename Type>
class Vector {
    using IndexType = int64_t;

   public:
    enum Orientation { Vertical = false, Horizontal = true };

    Vector() = default;

    explicit Vector(size_t length, Orientation orientation = Vertical)
        : data_(length, Type(0)), orientation_(orientation) {
        assert(length > 0);
    }

    Vector(const Type& element, size_t length, Orientation orientation = Vertical)
        : data_(length, element), orientation_(orientation) {
        assert(length > 0);
    }

    Vector(std::initializer_list<Type> list, Orientation orientation = Vertical)
        : data_(list), orientation_(orientation) {
        assert(list.size() > 0);
    }

    static Vector standart_basis(IndexType ind, size_t length, Orientation orientation = Vertical) {
        Vector result(length, orientation);
        result[ind] = Type(1.0);
        return result;
    }

    size_t size() const noexcept {
        return data_.size();
    }

    Orientation orientation() const noexcept {
        return orientation_;
    }

    Type& operator[](IndexType ind) {
        assert(ind >= 0 && ind < data_.size());
        return data_[ind];
    }

    const Type& operator[](IndexType ind) const {
        assert(ind >= 0 && ind < data_.size());
        return data_[ind];
    }

    Vector& operator+=(const Vector& rhs) noexcept {
        assert(data_.size() == rhs.size() && orientation_ == rhs.orientation_);
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] += rhs[ind];
        }
        return *this;
    }

    friend Vector operator+(const Vector& lhs, const Vector& rhs) noexcept {
        Vector result = lhs;
        result += rhs;
        return result;
    }

    Vector& operator-=(const Vector& rhs) noexcept {
        assert(data_.size() == rhs.size() && orientation_ == rhs.orientation_);
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] -= rhs[ind];
        }
        return *this;
    }

    friend Vector operator-(const Vector& lhs, const Vector& rhs) noexcept {
        Vector result = lhs;
        result -= rhs;
        return result;
    }

    Vector& operator*=(const Type& rhs) noexcept {
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] *= rhs;
        }
        return *this;
    }

    Vector& operator/=(const Type& rhs) noexcept {
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] /= rhs;
        }
        return *this;
    }

    Vector operator-() const noexcept {
        Vector result = *this;
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            result[ind] = -result[ind];
        }
        return result;
    }

    friend bool operator==(const Vector& lhs, const Vector& rhs) noexcept {
        return lhs.data_ == rhs.data_ && lhs.orientation_ == rhs.orientation_;
    }

    friend bool operator!=(const Vector& lhs, const Vector& rhs) noexcept {
        return !(lhs == rhs);
    }

    Vector& transpose() noexcept {
        if (orientation_ == Vertical) {
            orientation_ = Horizontal;
        } else {
            orientation_ = Vertical;
        }
        return *this;
    }

    friend std::istream& operator>>(std::istream& in, Vector<Type>& v) noexcept {
        size_t sz;
        in >> sz;
        v = Vector<Type>(sz);
        for (size_t ind = 0; ind < sz; ++ind) {
            in >> v[ind];
        }
        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const Vector<Type>& v) noexcept {
        char delimetr = (v.orientation_ == Vertical) ? '\n' : ' ';
        for (size_t ind = 0; ind < v.size(); ++ind) {
            out << v[ind];
            if (ind + 1 < v.size()) {
                out << delimetr;
            }
        }
        return out;
    }

   private:
    std::vector<Type> data_;
    Orientation orientation_;
};

template <typename Type>
Vector<Type> operator/(const Vector<Type>& lhs, const Type& rhs) noexcept {
    Vector<Type> result = lhs;
    result /= rhs;
    return result;
}

template <typename Type>
Vector<Type> operator*(const Vector<Type>& lhs, const Type& rhs) noexcept {
    Vector<Type> result = lhs;
    result *= rhs;
    return result;
}

template <typename Type>
Vector<Type> operator*(const Type& lhs, const Vector<Type>& rhs) noexcept {
    Vector<Type> result = rhs;
    result *= lhs;
    return result;
}

template <typename Type>
Vector<Type> transpose(const Vector<Type>& v) {
    Vector<Type> result = v;
    return result.transpose();
}

template <typename Type>
Type dot_product(const Vector<Type>& lhs, const Vector<Type>& rhs) noexcept;

template <>
inline long double dot_product(const Vector<long double>& lhs, const Vector<long double>& rhs) noexcept {
    assert(lhs.size() == rhs.size() && lhs.orientation() == Vector<long double>::Orientation::Horizontal &&
           rhs.orientation() == Vector<long double>::Orientation::Vertical);
    long double res = 0;
    for (size_t ind = 0; ind < lhs.size(); ++ind) {
        res += lhs[ind] * rhs[ind];
    }
    return res;
}

template <>
inline Complex dot_product(const Vector<Complex>& lhs, const Vector<Complex>& rhs) noexcept {
    assert(lhs.size() == rhs.size() && lhs.orientation() == Vector<Complex>::Orientation::Horizontal &&
           rhs.orientation() == Vector<Complex>::Orientation::Vertical);
    Complex res = 0;
    for (size_t ind = 0; ind < lhs.size(); ++ind) {
        res += conjugate(lhs[ind]) * rhs[ind];
    }
    return res;
}

template <typename Type>
long double abs(Vector<Type> v) noexcept;

template <>
inline long double abs(Vector<long double> v) noexcept {
    if (v.orientation() == Vector<long double>::Orientation::Vertical) {
        long double res = dot_product(transpose(v), v);
        return sqrtl(res);
    }
    long double res = dot_product(v, transpose(v));
    return sqrtl(res);
}

template <>
inline long double abs(Vector<Complex> v) noexcept {
    if (v.orientation() == Vector<Complex>::Orientation::Vertical) {
        long double res = dot_product(transpose(v), v).Re();
        return sqrtl(res);
    }
    long double res = dot_product(v, transpose(v)).Re();
    return sqrtl(res);
}

}  // namespace svd_computation
