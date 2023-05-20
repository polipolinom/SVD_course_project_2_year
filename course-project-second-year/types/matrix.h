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
        assert(data.size() > 0 && data.begin()[0].size() > 0);
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
        assert(height > 0 && width > 0);
        data_ = std::vector<Type>(height * width, Type(0));
    }

    size_t height() const noexcept {
        return height_;
    }

    size_t width() const noexcept {
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
        assert(row < height_ && row >= 0 && column >= 0 && column < width());
        return data_[row * width() + column];
    }

    const Type& operator()(IndexType row, IndexType column) const noexcept {
        assert(row < height_ && row >= 0 && column >= 0 && column <= width());
        return data_[row * width() + column];
    }

    static Matrix identity(size_t n) {
        assert(n > 0);
        Matrix result(n, n);
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = Type(1);
        }
        return result;
    }

    static Matrix diagonal(std::initializer_list<Type> diagonal) noexcept {
        assert(diagonal.size() > 0);

        Matrix result(diagonal.size(), diagonal.size());
        for (size_t ind = 0; ind < diagonal.size(); ++ind) {
            result(ind, ind) = diagonal.begin()[ind];
        }
        return result;
    }

    static Matrix from_vectors(std::vector<Vector<Type>> v) noexcept {
        assert(!v.empty());

        if (v[0].orientation() == Vector<Type>::Orientation::Horizontal) {
            Matrix result = Matrix(v.size(), v[0].size());
            for (size_t i = 0; i < v.size(); ++i) {
                assert(v[0].size() == v[i].size() && v[i].orientation() == v[0].orientation());
                for (size_t j = 0; j < v[i].size(); ++j) {
                    result(i, j) = v[i][j];
                }
            }
            return result;
        }

        Matrix result = Matrix(v[0].size(), v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            assert(v[0].size() == v[i].size() && v[i].orientation() == v[0].orientation());
            for (size_t j = 0; j < v[i].size(); ++j) {
                result(j, i) = v[i][j];
            }
        }
        return result;
    }

    Matrix& operator+=(const Matrix& rhs) noexcept {
        assert(rhs.height_ == height_ && rhs.width() == width());
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
        assert(rhs.height_ == height_ && rhs.width() == width());
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
        assert(rhs.height_ == width());

        Matrix result(height_, rhs.width());
        for (size_t row = 0; row < height_; ++row) {
            for (size_t column = 0; column < rhs.width(); ++column) {
                for (size_t k = 0; k < rhs.height_; ++k) {
                    result(row, column) += (*this)(row, k) * rhs(k, column);
                }
            }
        }
        *this = std::move(result);

        return *this;
    }

    friend Matrix operator*(const Matrix& lhs, const Matrix& rhs) noexcept {
        Matrix result = lhs;
        result *= rhs;
        return result;
    }

    friend Vector<Type> operator*(const Matrix& lhs, const Vector<Type>& rhs) noexcept {
        assert(rhs.size() == lhs.width() && rhs.orientation() == Vector<Type>::Orientation::Vertical);

        Vector<Type> result(rhs.size());
        for (size_t row = 0; row < lhs.height(); ++row) {
            for (size_t column = 0; column < rhs.size(); ++column) {
                result[row] += lhs(row, column) * rhs[column];
            }
        }

        return result;
    }

    friend Vector<Type> operator*(const Vector<Type>& lhs, const Matrix& rhs) noexcept {
        assert(lhs.size() == rhs.height() && lhs.orientation() == Vector<Type>::Orientation::Horizontal);

        Vector<Type> result(Type(0.0), rhs.size(), Vector<Type>::Orientation::Horizontal);

        for (size_t column = 0; column < rhs.width(); ++column) {
            for (size_t row = 0; row < lhs.size(); ++row) {
                result[column] += rhs[row] * rhs(row, column);
            }
        }

        return result;
    }

    Matrix& operator*=(const Type& rhs) noexcept {
        for (size_t row = 0; row < height_; ++row) {
            for (size_t column = 0; column < width(); ++column) {
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
        for (size_t row = 0; row < height_; ++row) {
            for (size_t column = 0; column < width(); ++column) {
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
            for (size_t column = 0; column < width(); column++) {
                (*this)(row, column) = -(*this)(row, column);
            }
        }
        return *this;
    }

    Matrix& transpose() noexcept {
        Matrix result(width(), height_);
        for (size_t row = 0; row < width(); ++row) {
            for (size_t column = 0; column < height_; ++column) {
                result(row, column) = (*this)(column, row);
            }
        }

        *this = std::move(result);
        return *this;
    }

    Matrix& conjugate() noexcept;

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
                ss.flags(out.flags());
                ss.precision(out.precision());
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

template <typename Type>
Matrix<Type> transpose(const Matrix<Type>& A) noexcept {
    Matrix<Type> result(A.width(), A.height());
    for (size_t row = 0; row < A.width(); ++row) {
        for (size_t column = 0; column < A.height(); ++column) {
            result(row, column) = A(column, row);
        }
    }
    return result;
}

template <>
inline Matrix<Complex>& Matrix<Complex>::conjugate() noexcept {
    transpose();
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i].conjugate();
    }

    return *this;
}

template <>
inline Matrix<long double>& Matrix<long double>::conjugate() noexcept {
    return transpose();
}

template <typename Type>
Matrix<Type> conjugate(const Matrix<Type>& A) noexcept;

template <>
inline Matrix<long double> conjugate(const Matrix<long double>& A) noexcept {
    return transpose(A);
}

template <>
inline Matrix<Complex> conjugate(const Matrix<Complex>& A) noexcept {
    Matrix<Complex> result = transpose(A);
    for (size_t i = 0; i < A.width(); ++i) {
        for (size_t j = 0; j < A.height(); ++j) {
            result(i, j).conjugate();
        }
    }

    return result;
}
}  // namespace svd_computation