#include "matrix.h"
#include <bits/stdc++.h>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

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
        phi_k = atan(2 * A.data[a][b] / (A.data[a][a] - A.data[b][b])) / 2.0;
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
  cout << "no. iterations: " << iter << "\n";
}

int main() {
  ifstream input_file("tests/1-4.json");
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
  long double e = data["eps"].get<long double>();

  Matrix MA = Matrix(A);
  Jacobi(MA, e);
}
