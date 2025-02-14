#include "../lib/linalg.h"
#include "../lib/matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(Matrix &A) {
  if (A.rows != A.cols) {
    throw std::invalid_argument("Matrix must be square");
  }
  std::size_t n = A.rows;
  Matrix L(n, n);
  Matrix U(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = i; k < n; ++k) {
      long double sum = 0.0;
      for (std::size_t j = 0; j < i; ++j) {
        sum += L.data[i][j] * U.data[j][k];
      }
      U.data[i][k] = A.data[i][k] - sum;
    }
    for (std::size_t k = i; k < n; ++k) {
      if (i == k) {
        L.data[i][i] = 1.0;
      } else {
        long double sum = 0.0;
        for (std::size_t j = 0; j < i; ++j) {
          sum += L.data[k][j] * U.data[j][i];
        }
        L.data[k][i] = (A.data[k][i] - sum) / U.data[i][i];
      }
    }
  }
  return std::make_pair(L, U);
}

std::vector<long double> forwardSubstitution(const Matrix &L,
                                             const std::vector<long double> &b) {
  std::vector<long double> x(L.rows, 0.0);
  for (std::size_t i = 0; i < L.rows; ++i) {
    long double sum = 0.0;
    for (std::size_t j = 0; j < i; ++j) {
      sum += L.data[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / L.data[i][i];
  }
  return x;
}

std::vector<long double> backwardSubstitution(const Matrix &U,
                                         const std::vector<long double> &y) {
  std::vector<long double> x(U.rows, 0.0);
  for (std::size_t i = U.rows - 1; i != std::string::npos; --i) {
    long double sum = 0.0;
    for (std::size_t j = i + 1; j < U.rows; ++j) {
      sum += U.data[i][j] * x[j];
    }
    x[i] = (y[i] - sum) / U.data[i][i];
  }
  return x;
}

Matrix solveLS(Matrix &A, const std::vector<long double> &b) {
    if (A.rows != A.cols || A.rows != b.size()) {
        throw std::invalid_argument("Matrix must be square and have the same number of rows as the vector");
    }
    auto [L, U] = LUDecomposition(A);
    std::vector<long double> y = forwardSubstitution(L, b);
    std::vector<long double> x = backwardSubstitution(U, y);
    return Matrix(x);
}
