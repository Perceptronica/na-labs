#include <bits/stdc++.h>
#include <cmath>
#include <nlohmann/json.hpp>
#include <regex>

using namespace std;
using json = nlohmann::json;

double f(double x) { return exp(x); }

string format_polynom(const string &polynom_str) {
  string result = polynom_str;
  result = regex_replace(result, regex("--"), "+");
  result = regex_replace(result, regex("\\+-"), "-");
  result = regex_replace(result, regex("\\+"), " + ");
  result = regex_replace(result, regex("\\* "), "*");
  result = regex_replace(result, regex("\\s+"), " ");
  result = regex_replace(result, regex("^\\s+|\\s+$"), "");
  return result;
}

pair<string, double> Lagrangia(const vector<double> &x, const vector<double> &y,
                               const pair<double, double> &test_point) {
  assert(x.size() == y.size());

  string polynom_str = "L(x) =";
  double polynom_test_value = 0;

  for (size_t i = 0; i < x.size(); ++i) {
    string cur_enum_str;
    double cur_enum_test = 1;
    double cur_denom = 1;

    for (size_t j = 0; j < x.size(); ++j) {
      if (i == j)
        continue;

      cur_enum_str += "(x-" + to_string(x[j]) + ")";
      cur_enum_test *= (test_point.first - x[j]);
      cur_denom *= (x[i] - x[j]);
    }

    polynom_str += "+" + to_string(y[i] / cur_denom) + cur_enum_str;
    polynom_test_value += y[i] * cur_enum_test / cur_denom;
  }

  return {format_polynom(polynom_str),
          abs(polynom_test_value - test_point.second)};
}

pair<string, double> Newton(const vector<double> &x, const vector<double> &y,
                            const pair<double, double> &test_point) {
  assert(x.size() == y.size());

  size_t n = x.size();
  vector<double> coefs(y);

  for (size_t i = 1; i < n; ++i) {
    for (size_t j = n - 1; j >= i; --j) {
      coefs[j] = (coefs[j] - coefs[j - 1]) / (x[j] - x[j - i]);
    }
  }

  string polynom_str = "P(x) = ";
  double polynom_test_value = 0;
  string cur_multipliers_str;
  double cur_multipliers = 1;

  for (size_t i = 0; i < n; ++i) {
    polynom_test_value += cur_multipliers * coefs[i];

    if (i == 0) {
      polynom_str += to_string(coefs[i]);
    } else {
      polynom_str += "+" + cur_multipliers_str + "*" + to_string(coefs[i]);
    }

    cur_multipliers *= (test_point.first - x[i]);
    cur_multipliers_str += "(x-" + to_string(x[i]) + ")";
  }

  return {format_polynom(polynom_str),
          abs(polynom_test_value - test_point.second)};
}

int main() {
  ifstream input_file("tests/3-1.json");
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

  vector<double> x_a = data["x_a"].get<vector<double>>();
  vector<double> x_b = data["x_b"].get<vector<double>>();
  double x_test = data["x*"].get<double>();

  vector<double> y_a(x_a.size());
  vector<double> y_b(x_b.size());
  for (size_t i = 0; i < x_a.size(); ++i) {
    y_a[i] = f(x_a[i]);
  }
  for (size_t i = 0; i < x_b.size(); ++i) {
    y_b[i] = f(x_b[i]);
  }

  pair<double, double> test_point = {x_test, f(x_test)};

  auto [lagr_poly_a, lagr_err_a] = Lagrangia(x_a, y_a, test_point);
  auto [lagr_poly_b, lagr_err_b] = Lagrangia(x_b, y_b, test_point);
  auto [newton_poly_a, newton_err_a] = Newton(x_a, y_a, test_point);
  auto [newton_poly_b, newton_err_b] = Newton(x_b, y_b, test_point);

  cout << "Для набора x_a:\n";
  cout << lagr_poly_a << "\nΔL(" << x_test << ") = " << lagr_err_a << endl;
  cout << newton_poly_a << "\nΔP(" << x_test << ") = " << newton_err_a << endl
       << endl;

  cout << "Для набора x_b:\n";
  cout << "Лагранж: " << lagr_poly_b << "\nΔL(" << x_test << ") = " << lagr_err_b << endl;
  cout << "Ньютон: " << newton_poly_b << "\nΔP(" << x_test << ") = " << newton_err_b << endl;
}
