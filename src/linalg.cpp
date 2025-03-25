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

std::vector<long double> createDiags(const Matrix &M, std::size_t x,
                                     std::size_t y, std::size_t n) {
  std::vector<long double> diag(n);
  for (std::size_t i = 0; i < n; ++i) {
    diag[i] = M.data[x][y];
    x++;
    y++;
  }
  return diag;
}

Matrix solveTDMA(const Matrix &A, const Matrix &B) {
  if (A.rows != A.cols || A.rows != B.rows) {
    throw std::invalid_argument(
        "Matrix must be square and have the same number of rows as the vector");
  }
  std::vector<long double> a = createDiags(A, 1, 0, A.rows - 1);
  std::vector<long double> b = createDiags(A, 0, 0, A.rows);
  std::vector<long double> c = createDiags(A, 0, 1, A.rows - 1);
  std::size_t n = A.rows;
  std::vector<long double> x(n);
  std::vector<long double> cp(n), dp(n);
  a.insert(a.begin(), 0);
  c.push_back(0);
  cp[0] = c[0] / b[0];
  dp[0] = B.data[0][0] / b[0];

  for (size_t i = 1; i < n; ++i) {
    double denominator = b[i] - a[i] * cp[i - 1];
    if (denominator == 0.0) {
      throw std::runtime_error("Division by zero in forward pass!");
    }
    if (i < n - 1) {
      cp[i] = c[i] / denominator;
    }
    dp[i] = (B.data[i][0] - a[i] * dp[i - 1]) / denominator;
  }

  x[n - 1] = dp[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    x[i] = dp[i] - cp[i] * x[i + 1];
  }
  Matrix ans = Matrix(x);
  return ans;
}

std::pair<Matrix, int> simpleIterationImpl(const Matrix &A, const Matrix &B,
                                           long double epsilon,
                                           int maxIterations, long double tau) {
  std::size_t n = A.rows;
  std::vector<long double> x(n, 0.0);
  std::vector<long double> x_prev(n);
  std::size_t i = 0;
  for (; i < maxIterations; ++i) {
    x_prev = x;
    for (std::size_t j = 0; j < n; ++j) {
      double Ax = 0.0;
      for (std::size_t k = 0; k < n; ++k) {
        Ax += A.data[j][k] * x_prev[k];
      }
      x[j] = x_prev[j] - tau * Ax + tau * B.data[j][0];
    }
    double diff = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      diff += (x[j] - x_prev[j]) * (x[j] - x_prev[j]);
    }
    diff = sqrt(diff);
    if (diff < epsilon) {
      break;
    }
  }
  return std::make_pair(Matrix(x), i);
}

long double vectorDiffNorm(std::vector<long double>& x, std::vector<long double>& y) {
    long double norm = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        norm += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return std::sqrt(norm);
}

std::pair<Matrix, int> Seidel(const Matrix &A, const Matrix &b, long double eps, int maxIters) {
  std::size_t n = A.rows;
  std::vector<long double> x(n, 0.0);
  std::vector<long double> x_prev(n);
  int iter = 0;
  for (; iter < maxIters; ++iter) {
    x_prev = x;
    for (size_t i = 0; i < n; ++i) {
      double sum1 = 0.0;
      for (size_t j = 0; j < i; ++j) {
        sum1 += A.data[i][j] * x[j];
      }
      double sum2 = 0.0;
      for (size_t j = i + 1; j < n; ++j) {
        sum2 += A.data[i][j] * x_prev[j];
      }
      x[i] = (b.data[i][0] - sum1 - sum2) / A.data[i][i];
    }

    if (vectorDiffNorm(x, x_prev) < eps) {
      break;
    }
  }
  return std::make_pair(Matrix(x), iter);
}
