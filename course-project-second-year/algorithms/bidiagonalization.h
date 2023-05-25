#pragma once

#include <cmath>
#include <vector>

#include "../types/matrix.h"
#include "../types/vector.h"
#include "../utils/complement_orthobase.h"
#include "../utils/set_values.h"

namespace svd_computation {
template <typename Type>
Matrix<long double> bidiagonalize(const Matrix<Type>& A, Matrix<Type>* left_basis = nullptr,
                                  Matrix<Type>* right_basis = nullptr,
                                  const long double eps = constants::DEFAULT_EPSILON) {
    Matrix<Type> B = A;
    Matrix<long double> result(A.height(), A.width());
    if (left_basis != nullptr) {
        *left_basis = Matrix<Type>::identity(A.height());
    }
    if (right_basis != nullptr) {
        *right_basis = Matrix<Type>::identity(A.width());
    }

    for (size_t ind = 0; ind < std::min(A.height(), A.width()); ++ind) {
        left_reflection(B, ind, ind, left_basis, eps);
        if (ind + 1 < A.width()) {
            right_reflection(B, ind, ind + 1, right_basis, eps);
        }
    }

    if (left_basis != nullptr) {
        (*left_basis).conjugate();
    }

    return B;
}

namespace details {
Matrix<long double> bidiagonalize_with_right_basis(Matrix<long double>& A, const Matrix<long double>& right_basis,
                                                   const long double eps) {
    using Matrix = Matrix<long double>;

    Matrix basis = Matrix::identity(A.height());

    Matrix A1 = A * right_basis;

    std::vector<Vector<long double>> U;

    long double alpha = 0;
    long double beta = 0;
    Vector<long double> u(0, A.height());
    for (size_t ind = 0; ind < A.height(); ++ind) {
        alpha = abs(A * right_basis.column(ind) - beta * u);
        if (alpha <= eps) {
            u = add_one_vector(U, A.height(), Vector<long double>::Orientation::Vertical, eps);
        } else {
            u = (A1.column(ind) - beta * u) / alpha;
        }
        for (size_t i = 0; i < u.size(); ++i) {
            basis(i, ind) = u[i];
        }

        U.emplace_back(u);

        if (ind + 1 < A.width()) {
            beta = dot_product(transpose(u), A1.column(ind + 1));
        }
    }
    return basis;
}
}  // namespace details
}  // namespace svd_computation
