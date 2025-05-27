#include "matrix.h"
#include "norms.h"
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

void simple_iters(const Matrix &A, const Matrix &B,
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

int main() {
  ifstream input_file("tests/1-3.json");
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
  long double e = data["eps"].get<long double>();

  Matrix MA = Matrix(A);
  Matrix mb = Matrix(b);
  cout << "Simple iteration method:\n";
  simple_iters(MA, mb, e);
  cout << "\nSeidel's method:\n";
  Seidel(MA, mb, e);
  cout << "\n";
}
