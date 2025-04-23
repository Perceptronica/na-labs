#include "../lib/linalg.h"
#include "../lib/matrix.h"
#include <climits>
#include <cmath>

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
  std::cout << L << "\n";
  std::cout << U << "\n";
  long double det = 1.0;
  for (std::size_t i = 0; i < U.rows; ++i) {
    det *= U.data[i][i];
  }
  std::cout << "'smart' determinant: " << det << "\n";
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
  std::vector<long double> P(n), Q(n);
  a.insert(a.begin(), 0);
  c.push_back(0);
  P[0] = (-1) * c[0] / b[0];
  Q[0] = B.data[0][0] / b[0];
  for (std::size_t i = 1; i < n; ++i) {
    long double denom = b[i] + (a[i] * P[i - 1]);
    P[i] = ((-1) * c[i]) / denom;
    Q[i] = (B.data[i][0] - (a[i] * Q[i - 1])) / denom;
  }
  x[n - 1] = Q[n - 1];
  for (std::size_t i = n - 2; i > 0; --i) {
    x[i] = P[i] * x[i + 1] + Q[i];
  }
  x[0] = P[0] * x[1] + Q[0];
  Matrix ans = Matrix(x);
  return ans;
}

void simpleIterationImpl(const Matrix &A, const Matrix &B,
                         long double epsilon) {
  Matrix alpha = A;
  Matrix beta = B;
  for (std::size_t i = 0; i < A.rows; ++i) {
    long double c = alpha.data[i][i];
    for (std::size_t j = 0; j < A.cols; ++j) {
      alpha.data[i][j] /= c;
      alpha.data[i][j] *= -1;
    }
    beta.data[i][0] /= c;
    alpha.data[i][i] = 0;
  }
  // std::cout << norm(3, alpha) << "\n";
  Matrix x = beta;
  Matrix x_next = x;
  int i = 1;
  do {
    i++;
    x = x_next;
    x_next = beta + (alpha * x);
  } while (norm(3, x_next - x) > epsilon);
  std::cout << x_next << "\n";
  std::cout << "no. iterations: " << i << "\n";
}

long double norm(std::size_t mode, const Matrix &V) {
  long double norm = 0.0;
  if (V.cols == 1) {
    if (mode == 1) {
      for (std::size_t i = 0; i < V.rows; ++i) {
        norm += fabs(V.data[i][0]);
      }
    }
    if (mode == 2) {
      for (std::size_t i = 0; i < V.rows; ++i) {
        norm += (V.data[i][0] * V.data[i][0]);
      }
      norm = sqrt(norm);
    }
    if (mode == 3) {
      norm = INT_MIN;
      for (std::size_t i = 0; i < V.rows; ++i) {
        norm = fmax(norm, fabs(V.data[i][0]));
      }
    }
  } else {
    if (mode == 1) {
      for (std::size_t j = 0; j < V.cols; ++j) {
        long double c = INT_MIN;
        for (std::size_t i = 0; i < V.rows; ++i) {
          c += fabs(V.data[i][j]);
        }
        norm = fmax(norm, c);
      }
    }
    if (mode == 2) {
      for (std::size_t j = 0; j > V.cols; ++j) {
        for (std::size_t i = 0; i < V.rows; ++i) {
          norm += (V.data[i][j] * V.data[i][j]);
        }
      }
      norm = sqrt(norm);
    }
    if (mode == 3) {
      norm = INT_MIN;
      for (std::size_t i = 0; i < V.rows; ++i) {
        long double c = 0;
        for (std::size_t j = 0; j < V.cols; ++j) {
          c += fabs(V.data[i][j]);
        }
        norm = fmax(norm, c);
      }
    }
  }
  return norm;
}

void Seidel(const Matrix &a, const Matrix &b, long double epsilon) {
  Matrix alpha = a;
  Matrix beta = b;
  for (std::size_t i = 0; i < a.rows; ++i) {
    long double c = alpha.data[i][i];
    for (std::size_t j = 0; j < a.cols; ++j) {
      alpha.data[i][j] /= c;
      alpha.data[i][j] *= -1;
    }
    beta.data[i][0] /= c;
    alpha.data[i][i] = 0;
  }
  // std::cout << norm(3, alpha) << "\n";
  Matrix B(alpha.rows, alpha.cols);
  Matrix C(alpha.rows, alpha.cols);
  for (std::size_t i = 0; i < B.rows; ++i) {
    for (std::size_t j = 0; j < B.rows; ++j) {
      if (j < i) {
        B.data[i][j] = alpha.data[i][j];
      } else {
        C.data[i][j] = alpha.data[i][j];
      }
    }
  }
  Matrix x = beta;
  Matrix x_next = x;
  Matrix E(B.rows, B.cols);
  for (std::size_t i = 0; i < E.rows; ++i) {
    E.data[i][i] = 1;
  }
  Matrix D(B.rows, B.cols);
  int i = 1;
  do {
    i++;
    x = x_next;
    D = (E - B).inverse();
    x_next = D * C * x + D * beta;
  } while (norm(3, x_next - x) > epsilon);
  std::cout << x_next << "\n";
  std::cout << "no. iterations: " << i << "\n";
}

void Jacobi(const Matrix &Q, long double epsilon) {
  Matrix A = Q;
  std::vector<Matrix> U_vec;
  long double err = epsilon;
  int iter = 0;
  do {
    ++iter;
    err = 0.0;
    long double m = 0.0;
    std::size_t a = 1;
    std::size_t b = 1;
    for (std::size_t i = 0; i < A.rows; ++i) {
      for (std::size_t j = i + 1; j < A.cols; ++j) {
        if (fabs(A.data[i][j]) > m) {
          m = fabs(A.data[i][j]);
          a = i;
          b = j;
        }
      }
    }
    if (m < epsilon) {
        break;
    }
    long double phi_k;
    if (A.data[a][a] == A.data[b][b]) {
        phi_k = M_PI / 4;
    } else {
        phi_k = atan((2 * m) / (A.data[a][a] - A.data[b][b])) / 2.0;
    }
    Matrix U(A.rows, A.cols);
    for (std::size_t i = 0; i < U.rows; ++i) {
      U.data[i][i] = 1;
    }
    U.data[a][b] = -sin(phi_k);
    U.data[b][a] = sin(phi_k);
    U.data[a][a] = cos(phi_k);
    U.data[b][b] = cos(phi_k);
    // std::cout << U << "\n";
    U_vec.push_back(U);
    U.transpose();
    A = U * A;
    U.transpose();
    A = A * U;
    for (std::size_t i = 0; i < A.rows; ++i) {
      for (std::size_t j = i + 1; j < A.cols; ++j) {
        err += (A.data[i][j] * A.data[i][j]);
      }
    }
  } while (sqrtl(err) > epsilon);
  std::cout << "Eigenvalues:\n";
  for (std::size_t i = 0; i < A.rows; ++i) {
    std::cout << A.data[i][i] << "\n";
  }
  Matrix res = U_vec[0];
  for (std::size_t i = 1; i < U_vec.size(); ++i) {
    res = res * U_vec[i];
  }
  std::cout << "Eigenvectors:\n";
  for (std::size_t i = 0; i < iter - 1; ++i) {
    std::cout << "U" << i << (i == iter - 2 ? "" : "*");
  }
  std::cout << "=\n" << res << "\n";
  for (std::size_t i = 0; i < res.cols; ++i) {
    std::cout << i + 1 << ":\n";
    for (std::size_t j = 0; j < res.rows; ++j) {
      std::cout << res.data[j][i] << "\n";
    }
  }
}
