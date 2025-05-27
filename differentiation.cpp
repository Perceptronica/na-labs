#include <bits/stdc++.h>
#include <cassert>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

double df(const vector<double> &x, const vector<double> &y, double x_) {
  assert(x.size() == y.size());
  size_t i = 0;
  for (size_t interval = 0; interval < x.size() - 1; ++interval) {
    if (x[interval] <= x_ && x_ < x[interval + 1]) {
      i = interval;
      break;
    }
  }
  double a1 = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
  double a2 = ((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - a1) /
              (x[i + 2] - x[i]) * (2 * x_ - x[i] - x[i + 1]);
  return a1 + a2;
}

double d2f(const vector<double> &x, const vector<double> &y, double x_) {
  assert(x.size() == y.size());
  size_t i = 0;
  for (size_t interval = 0; interval < x.size() - 1; ++interval) {
    if (x[interval] <= x_ && x_ < x[interval + 1]) {
      i = interval;
      break;
    }
  }
  double num = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) -
               (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
  return 2 * num / (x[i + 2] - x[i]);
}

int main() {
  ifstream input_file("tests/3-4.json");
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

  vector<double> x = data["x"].get<vector<double>>();
  vector<double> y = data["y"].get<vector<double>>();
  double x_ = data["x*"].get<double>();

  cout << "f'(" << x_ << ") = " << df(x, y, x_) << endl;
  cout << "f''(" << x_ << ") = " << d2f(x, y, x_) << endl;
}
