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
    Vector() = default;

    explicit Vector(size_t length, bool horizontal = false) : data_(length, Type(0)), is_horizontal_(horizontal) {}

    Vector(const Type& element, size_t length, bool horizontal = false)
        : data_(length, Type(element)), is_horizontal_(horizontal) {}

    Vector(std::initializer_list<Type> list, bool horizontal = false) : data_(list), is_horizontal_(horizontal) {}

    size_t size() const noexcept {
        return data_.size();
    }

    bool is_vertical() const noexcept {
        return !is_horizontal_;
    }

    bool is_horizontal() const noexcept {
        return is_horizontal_;
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
        assert(data_.size() == rhs.size());
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
        assert(data_.size() == rhs.size());
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

    friend Vector operator*(const Type& lhs, const Vector& rhs) noexcept {
        Vector result = rhs;
        result *= lhs;
        return result;
    }

    Vector operator/=(const Type& rhs) noexcept {
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] /= rhs;
        }
        return *this;
    }

    friend Vector operator/(const Vector& lhs, const Type& rhs) noexcept {
        Vector result = lhs;
        result /= rhs;
        return result;
    }

    Vector operator-() const noexcept {
        Vector result = *this;
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            result[ind] = -result[ind];
        }
        return result;
    }

    friend bool operator==(const Vector& lhs, const Vector& rhs) noexcept {
        assert(lhs.size() == rhs.size());
        for (size_t ind; ind < lhs.size(); ++ind) {
            if (lhs[ind] != rhs[ind]) {
                return false;
            }
        }
        return true;
    }

    friend bool operator!=(const Vector& lhs, const Vector& rhs) noexcept {
        return !(lhs == rhs);
    }

    Vector transpose() const noexcept {
        Vector res = *this;
        res.is_horizontal_ = true;
        return res;
    }

    bool empty() const noexcept {
        return data_.empty();
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
        char delimetr = v.is_vertical() ? '\n' : ' ';
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
    bool is_horizontal_;
};

template <typename Type>
Type dot_product(Vector<Type> lhs, Vector<Type> rhs) noexcept;

template <>
long double dot_product(Vector<long double> lhs, Vector<long double> rhs) noexcept {
    assert(lhs.size() == rhs.size());
    assert(lhs.is_horizontal());
    assert(rhs.is_vertical());
    long double res = 0;
    for (size_t ind = 0; ind < lhs.size(); ++ind) {
        res += lhs[ind] * rhs[ind];
    }
    return res;
}

template <>
Complex dot_product(Vector<Complex> lhs, Vector<Complex> rhs) noexcept {
    assert(lhs.size() == rhs.size());
    assert(lhs.is_horizontal());
    assert(rhs.is_vertical());
    Complex res = 0;
    for (size_t ind = 0; ind < lhs.size(); ++ind) {
        res += lhs[ind].conjugate() * rhs[ind];
    }
    return res;
}

template <typename Type>
long double abs(Vector<Type> v) noexcept;

template <>
long double abs(Vector<long double> v) noexcept {
    return sqrtl(dot_product<long double>(v.transpose(), v));
}

template <>
long double abs(Vector<Complex> v) noexcept {
    long double res = dot_product<Complex>(v.transpose(), v).Re();
    return sqrtl(res);
}

}  // namespace svd_computation