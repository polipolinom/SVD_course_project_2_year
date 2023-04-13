#pragma once

#include <cassert>
#include <iostream>

namespace svd_computation {
template <typename Type>
class Vector {
    using IndexType = int64_t;

   public:
    Vector() = default;

    explicit Vector(size_t length) : data_(length, Type(0)), is_horizontal_(false) {}

    Vector(const Type& element, size_t length) : data_(length, Type(element)), is_horizontal_(false) {}

    Vector(std::initializer_list<Type> list) : data_(list), is_horizontal_(false) {}

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
        is_horizontal_ = true;
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
}  // namespace svd_computation
