#include "../lib/linalg.h"
#include "../lib/matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(const Matrix &A) {
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

Matrix forwardSubstitution(const Matrix &L, const Matrix &b) {
  std::vector<long double> x(L.rows, 0.0);
  for (std::size_t i = 0; i < L.rows; ++i) {
    long double sum = 0.0;
    for (std::size_t j = 0; j < i; ++j) {
      sum += L.data[i][j] * x[j];
    }
    x[i] = (b.data[i][0] - sum) / L.data[i][i];
  }
  return Matrix(x);
}

Matrix backwardSubstitution(const Matrix &U, const Matrix &y) {
  std::vector<long double> x(U.rows, 0.0);
  for (std::size_t i = U.rows - 1; i != std::string::npos; --i) {
    long double sum = 0.0;
    for (std::size_t j = i + 1; j < U.rows; ++j) {
      sum += U.data[i][j] * x[j];
    }
    x[i] = (y.data[i][0] - sum) / U.data[i][i];
  }
  return Matrix(x);
}

Matrix solveLU(const Matrix &A, const Matrix &b) {
  if (A.rows != A.cols || A.rows != b.rows) {
    throw std::invalid_argument(
        "Matrix must be square and have the same number of rows as the vector");
  }
  auto [L, U] = LUDecomposition(A);
  Matrix y = forwardSubstitution(L, b);
  Matrix x = backwardSubstitution(U, y);
  return x;
}
