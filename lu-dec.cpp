#include "matrix.h"
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

int p = 0;

tuple<Matrix, Matrix, Matrix> LU_dec(const Matrix &A) {
  size_t n = A.rows;

  Matrix P(n, n), L(n, n);
  Matrix U = A;

  for (size_t i = 0; i < n; ++i) { // делаем единичными
    P.data[i][i] = 1.0;
    L.data[i][i] = 1.0;
  }

  for (size_t i = 0; i < n;
       ++i) { // находим строку с максимальным эл-ом в текущем столбце
    size_t max_row = i;
    double max_value = abs(U.data[i][i]);

    for (size_t k = i + 1; k < n; ++k) {
      if (abs(U.data[k][i]) > max_value) {
        max_value = abs(U.data[k][i]);
        max_row = k;
      }
    }

    if (max_row !=
        i) { // если нашли строку с большим элементом, делаем перестановку
      U.data[i].swap(U.data[max_row]);
      P.data[i].swap(P.data[max_row]);
      for (size_t j = 0; j < i; ++j) {
        swap(L.data[i][j], L.data[max_row][j]);
      }
      ++p;
    }

    for (size_t j = i + 1; j < n; ++j) { // Gaussian elimination
      L.data[j][i] = U.data[j][i] / U.data[i][i];
      for (size_t k = i; k < n; ++k) {
        U.data[j][k] -= L.data[j][i] * U.data[i][k];
      }
    }
  }
  return make_tuple(P, L, U);
}

Matrix LU_solve(const Matrix &A, Matrix &b) {
  auto [P, L, U] = LU_dec(A);
  if (L.rows != b.rows || b.cols != 1) {
    throw std::runtime_error("Incompatible matrix dimensions");
  }

  Matrix b_permuted = P * b;

  int n = L.rows;
  vector<long double> y(n, 0);

  // Ly = b_permuted
  for (int i = 0; i < n; ++i) {
    y[i] = b_permuted.data[i][0];
    for (int j = 0; j < i; ++j) {
      y[i] -= L.data[i][j] * y[j];
    }
    y[i] /= L.data[i][i];
  }

  // Ux = y
  vector<long double> x(n, 0);
  for (int i = n - 1; i >= 0; --i) {
    x[i] = y[i];
    for (int j = i + 1; j < n; ++j) {
      x[i] -= U.data[i][j] * x[j];
    }
    x[i] /= U.data[i][i];
  }

  vector<vector<long double>> result(n, vector<long double>(1));
  for (int i = 0; i < n; ++i) {
    result[i][0] = x[i];
  }
  return Matrix(result);
}

long double det(Matrix &A) {
  int p_curr = p;
  auto [P, L, U] = LU_dec(A);
  long double x = 1;
  for (size_t i = 0; i < A.rows; ++i) {
    x *= U.data[i][i];
  }
  x *= (p_curr % 2) ? -1 : 1;
  return x;
}

Matrix inverse_matrix(const Matrix& A) {
    if (A.rows != A.cols) {
        throw std::invalid_argument("matrix must be square for inversion");
    }
    size_t n = A.rows;
    Matrix inv_A(n, n); // создаем матрицу n x n для результата
    for (size_t i = 0; i < n; ++i) {
        // создаем i-й базисный вектор
        std::vector<std::vector<long double>> e_data(n, std::vector<long double>(1, 0.0));
        e_data[i][0] = 1.0;
        Matrix e(e_data);

        // A * x = e_i
        Matrix x = LU_solve(A, e);
        for (size_t j = 0; j < n; ++j) {
            inv_A.data[j][i] = x.data[j][0];
        }
    }
    return inv_A;
}

int main() {
  ifstream input_file("tests/1-1.json");
  if (!input_file.is_open()) {
    cerr << "json file opening failure" << endl;
    return 1;
  }
  json data;
  try {
    input_file >> data;
  } catch (const exception &e) {
    cerr << "json file reading error: " << e.what() << endl;
    return 1;
  }

  vector<vector<long double>> A;
  for (const auto &row : data["A"]) {
    vector<long double> rowVec;
    for (double val : row) {
      rowVec.push_back(val);
    }
    A.push_back(rowVec);
  }
  vector<long double> b = data["b"].get<vector<long double>>();

  Matrix MA = Matrix(A);
  Matrix mb = Matrix(b);
  auto [P, L, U] = LU_dec(Matrix(A));
  cout << "P:\n" << P << "\nL:\n" << L << "\nU:\n" << U << "\n";
  cout << "x:\n" << LU_solve(MA, mb) << "\n";
  cout << "L*U:\n" << L * U << "\nP*A:\n" << P * A << "\n";
  cout << "det(A) = " << det(MA) << "\n";
  cout << "A^(-1):\n" << inverse_matrix(MA) << "\n";
  cout << "A*A^(-1):\n" << MA * inverse_matrix(MA) << "\n";
}
