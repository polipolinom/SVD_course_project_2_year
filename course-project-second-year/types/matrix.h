#pragma once

#include <string.h>

#include <sstream>
#include <vector>

#include "vector_and_matrix.h"

namespace svd_computation {
template <typename Type>
class Matrix {
   public:
    Matrix() {}

    Matrix(Type lamda, size_t n)
        : _data(std::vector<std::vector<Type>>(n, std::vector<Type>(n, 0))), _height(n), _width(n) {
        for (size_t ind = 0; ind < n; ind++) {
            _data[ind][ind] = lamda;
        }
    }

    Matrix(std::vector<std::vector<Type>> data) : _data(data), _height(data.size()) {
        if (data.size() > 0) {
            _width = data[0].size();
        } else {
            _width = 0;
        }
    }

    Matrix(Vector<Type> v)
        : _data(std::vector<std::vector<Type>>(v.size(), std::vector<Type>(1, 0))), _height(v.size()), _width(1) {
        for (size_t ind = 0; ind < _height; ind++) {
            _data[ind][0] = v[ind];
        }
    }

    explicit Matrix(size_t height, size_t width) : _height(height), _width(width) {
        _data = std::vector<std::vector<Type>>(height, std::vector<Type>(width, Type()));
    }

    size_t get_height() const {
        return _height;
    }

    size_t get_width() const {
        return _width;
    }

    std::vector<Type>& operator[](size_t ind) {
        return _data[ind];
    }

    const std::vector<Type>& operator[](size_t ind) const {
        return (const std::vector<Type>&)_data[ind];
    }

    Matrix operator+(const Matrix& rhs) const {
        if (rhs._height != _height || rhs._width != _width) {
            throw std::invalid_argument("Incorrected sizes in matrix summation");
        }

        Matrix ans(_data);
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < _width; j++) {
                ans[i][j] += rhs[i][j];
            }
        }
        return ans;
    }

    Matrix& operator+=(const Matrix& rhs) {
        *this = *this + rhs;
        return *this;
    }

    Matrix operator-(const Matrix& rhs) const {
        if (rhs._height != _height || rhs._width != _width) {
            throw std::invalid_argument("Incorrected sizes in matrix subtraction");
        }

        Matrix ans(_data);
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < _width; j++) {
                ans[i][j] -= rhs[i][j];
            }
        }
        return ans;
    }

    Matrix& operator-=(const Matrix& rhs) {
        *this = *this - rhs;
        return *this;
    }

    Matrix operator*(const Matrix& rhs) const {
        if (rhs._height != _width) {
            throw std::invalid_argument("Incorrected sizes in matrix multiplication");
        }

        Matrix ans(_height, rhs._width);
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < rhs._width; j++) {
                ans[i][j] = 0;
                for (size_t k = 0; k < _width; k++) {
                    ans[i][j] += _data[i][k] * rhs[k][j];
                }
            }
        }
        return ans;
    }

    Matrix operator*(const Vector<Type>& rhs) {
        if (rhs.sisze() != _width) {
            throw std::invalid_argument("Incorrected sizes in matrix ans vector multiplication");
        }

        Vector<Type> ans(0, rhs.size());
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < _width; j++) {
                ans[i] += _data[i][j] * rhs[j];
            }
        }
        return ans;
    }

    Matrix& operator*=(const Matrix& rhs) {
        *this = *this * rhs;
        return *this;
    }

    Matrix& operator*=(const Type& rhs) {
        *this = rhs * (*this);
        return *this;
    }

    Matrix& operator/=(const Type& rhs) {
        *this = *this / rhs;
        return *this;
    }

    Matrix T() const {
        Matrix ans(_width, _height);
        for (size_t i = 0; i < _width; i++) {
            for (size_t j = 0; j < _height; j++) {
                ans[i][j] = _data[j][i];
            }
        }
        return ans;
    }

   private:
    std::vector<std::vector<Type>> _data;
    size_t _height;
    size_t _width;
};

template <typename Type>
Matrix<Type> operator*(const Type& lambda, const Matrix<Type>& A) {
    Matrix<Type> ans = A;
    for (size_t i = 0; i < ans.get_height(); i++) {
        for (size_t j = 0; j < ans.get_width(); j++) {
            ans[i][j] *= lambda;
        }
    }
    return ans;
}

template <typename Type>
Matrix<Type> operator*(const Matrix<Type>& A, const Type& lambda) {
    Matrix<Type> ans = A;
    for (size_t i = 0; i < ans.get_height(); i++) {
        for (size_t j = 0; j < ans.get_width(); j++) {
            ans[i][j] /= lambda;
        }
    }
    return ans;
}

template <typename Type>
std::istream& operator>>(std::istream& in, Matrix<Type>& A) {
    size_t height, width;
    Matrix<Type> res(height, width);
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            in >> res[i][j];
        }
    }
    A = res;
    return in;
}

template <typename Type>
std::ostream& operator<<(std::ostream& out, const Matrix<Type>& A) {
    size_t height = A.get_height();
    size_t width = A.get_width();
    std::vector<std::vector<std::string>> A_str(height, std::vector<std::string>(width, ""));
    size_t max_len = 0;
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            std::stringstream ss;
            ss << A[i][j];
            A_str[i][j] = ss.str();
            max_len = std::max(max_len, ss.str().size());
        }
    }
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            out << A_str[i][j];
            for (int _ = 0; _ + A_str[i][j].size() < max_len; _++) {
                out << " ";
            }
            if (j + 1 < width) {
                out << " ";
            }
        }
        if (i + 1 < height) {
            out << "\n";
        }
    }
    return out;
}

}  // namespace svd_computation