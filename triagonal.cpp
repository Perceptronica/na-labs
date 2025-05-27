#include "matrix.h"
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

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

int main() {
  ifstream input_file("tests/1-2.json");
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
  cout << "x:\n" << solveTDMA(MA, mb) << "\n";
}
