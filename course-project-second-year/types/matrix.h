#pragma once

#include <string.h>

#include <sstream>
#include <vector>

#include "vector.h"

namespace svd_computation {
template <typename Type>
class Matrix {
    using IndexType = int64_t;

   public:
    Matrix() = default;

    Matrix(std::initializer_list<std::initializer_list<Type>> data) : height_(data.size()) {
        if (data.size() == 0) {
            return;
        }
        size_t width = data.begin()[0].size();
        data_ = std::vector<Type>(data.size() * width, Type(0));
        for (size_t row = 0; row < data.size(); ++row) {
            assert(data.begin()[row].size() == width);
            for (size_t column = 0; column < width; ++column) {
                data_[row * width + column] = data.begin()[row].begin()[column];
            }
        }
    }

    Matrix(size_t height, size_t width) : height_(height) {
        data_ = std::vector<Type>(height * width, Type(0));
    }

    inline size_t height() const noexcept {
        return height_;
    }

    inline size_t width() const noexcept {
        return data_.size() / height_;
    }

    Vector<Type> column(IndexType ind) const {
        assert(ind >= 0);
        assert(ind < width());
        Vector<Type> res(height_);
        for (size_t k = 0; k < height_; ++k) {
            res[k] = (*this)(k, ind);
        }
        return res;
    }

    Vector<Type> row(IndexType ind) const {
        assert(ind >= 0);
        assert(ind < height_);
        Vector<Type> res(width());
        for (size_t k = 0; k < width(); ++k) {
            res[k] = (*this)(ind, k);
        }
        return res;
    }

    Type& operator()(IndexType row, IndexType column) noexcept {
        assert(row < height_ && row >= 0 && column >= 0 && column < (*this).width());
        return data_[row * (*this).width() + column];
    }

    const Type& operator()(IndexType row, IndexType column) const noexcept {
        assert(row < height_ && row >= 0 && column >= 0 && column <= (*this).width());
        return data_[row * (*this).width() + column];
    }

    static Matrix ones(size_t n) {
        Matrix result(n, n);
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = Type(1);
        }
        return result;
    }

    static Matrix diagonal(std::initializer_list<Type> diagonal) noexcept {
        Matrix result(diagonal.size(), diagonal.size());
        for (size_t ind = 0; ind < diagonal.size(); ++ind) {
            result(ind, ind) = diagonal.begin()[ind];
        }
        return result;
    }

    static Matrix from_vectors(std::vector<Vector<Type>> v) noexcept {
        if (v.empty()) {
            Matrix result = Matrix(0, 0);
            return result;
        }

        if (v[0].is_horizontal()) {
            Matrix result = Matrix(v.size(), v[0].size());
            for (size_t i = 0; i < v.size(); ++i) {
                assert(v[0].size() == v[i].size());
                assert(v[i].is_horizontal());
                for (size_t j = 0; j < v[i].size(); ++j) {
                    result(i, j) = v[i][j];
                }
            }
            return result;
        }

        Matrix result = Matrix(v[0].size(), v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            assert(v[0].size() == v[i].size());
            assert(v[i].is_vertical());
            for (size_t j = 0; j < v[i].size(); ++j) {
                result(j, i) = v[i][j];
            }
        }
        return result;
    }

    Matrix& operator+=(const Matrix& rhs) noexcept {
        assert(rhs.height_ == height_ && rhs.width() == (*this).width());
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] += rhs.data_[ind];
        }
        return *this;
    }

    friend Matrix operator+(const Matrix& lhs, const Matrix& rhs) noexcept {
        Matrix result = lhs;
        result += rhs;
        return result;
    }

    Matrix& operator-=(const Matrix& rhs) noexcept {
        assert(rhs.height_ == height_ && rhs.width() == (*this).width());
        for (size_t ind = 0; ind < data_.size(); ++ind) {
            data_[ind] -= rhs.data_[ind];
        }
        return *this;
    }

    friend Matrix operator-(const Matrix& lhs, const Matrix& rhs) noexcept {
        Matrix result = lhs;
        result -= rhs;
        return result;
    }

    Matrix& operator*=(const Matrix& rhs) noexcept {
        assert(rhs.height_ == (*this).width());

        Matrix result(height_, rhs.width());
        for (size_t row = 0; row < height_; ++row) {
            for (size_t column = 0; column < rhs.width(); ++column) {
                for (size_t k = 0; k < rhs.height_; k++) {
                    result(row, column) += (*this)(row, k) * rhs(k, column);
                }
            }
        }
        *this = result;

        return *this;
    }

    friend Matrix operator*(const Matrix& lhs, const Matrix& rhs) noexcept {
        Matrix result = lhs;
        result *= rhs;
        return result;
    }

    friend Vector<Type> operator*(const Matrix& lhs, const Vector<Type>& rhs) noexcept {
        assert(rhs.size() == lhs.width());
        assert(rhs.is_vertical());

        Vector<Type> result(rhs.size());
        for (size_t row = 0; row < lhs.height(); ++row) {
            for (size_t column = 0; column < rhs.size(); ++column) {
                result[row] += lhs(row, column) * rhs[column];
            }
        }
        return result;
    }

    friend Vector<Type> operator*(const Vector<Type>& lhs, const Matrix& rhs) noexcept {
        assert(lhs.size() == rhs.height());
        assert(lhs.is_horizontal());

        Vector<Type> result(0, rhs.size());
        result.transpose();
        for (size_t column = 0; column < rhs.width(); ++column) {
            for (size_t row = 0; row < lhs.size(); ++row) {
                result[column] += rhs[row] * rhs(row, column);
            }
        }
        return result;
    }

    Matrix& operator*=(const Type& rhs) noexcept {
        for (size_t row = 0; row < height_; ++row) {
            for (size_t column = 0; column < (*this).width(); ++column) {
                (*this)(row, column) *= rhs;
            }
        }
        return *this;
    }

    friend Matrix operator*(const Type& lhs, const Matrix& rhs) noexcept {
        Matrix result = rhs;
        result *= lhs;
        return result;
    }

    Matrix& operator/=(const Type& rhs) noexcept {
        for (size_t row = 0; row < height_; row++) {
            for (size_t column = 0; column < (*this).width(); column++) {
                (*this)(row, column) /= rhs;
            }
        }
        return *this;
    }

    friend Matrix operator/(const Matrix& lhs, const Type& rhs) noexcept {
        Matrix result = lhs;
        result /= rhs;
        return result;
    }

    Matrix& operator-() const {
        for (size_t row = 0; row < height_; row++) {
            for (size_t column = 0; column < (*this).width(); column++) {
                (*this)(row, column) = -(*this)(row, column);
            }
        }
        return *this;
    }

    bool empty() const noexcept {
        return data_.empty();
    }

    Matrix transpose() const noexcept {
        Matrix result((*this).width(), height_);
        for (size_t row = 0; row < (*this).width(); ++row) {
            for (size_t column = 0; column < height_; ++column) {
                result(row, column) = (*this)(column, row);
            }
        }
        return result;
    }

    Matrix conjugate() const noexcept;

    friend std::istream& operator>>(std::istream& in, Matrix<Type>& A) noexcept {
        size_t height, width;
        Matrix<Type> res(height, width);
        for (size_t i = 0; i < height; i++) {
            for (size_t j = 0; j < width; j++) {
                in >> res(i, j);
            }
        }
        A = res;
        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix<Type>& A) noexcept {
        size_t height = A.height();
        size_t width = A.width();
        std::vector<std::vector<std::string>> A_str(height, std::vector<std::string>(width, ""));
        size_t max_len = 0;
        for (size_t i = 0; i < height; i++) {
            for (size_t j = 0; j < width; j++) {
                std::stringstream ss;
                ss << A(i, j);
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

   private:
    std::vector<Type> data_;
    size_t height_ = 0;
};

template <>
Matrix<Complex> Matrix<Complex>::conjugate() const noexcept {
    Matrix result = transpose();
    for (size_t i = 0; i < width(); ++i) {
        for (size_t j = 0; j < height(); ++j) {
            result(i, j) = result(i, j).conjugate();
        }
    }
    return result;
}

template <>
Matrix<long double> Matrix<long double>::conjugate() const noexcept {
    return transpose();
}
}  // namespace svd_computation