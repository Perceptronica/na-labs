#include <bits/stdc++.h>
#include <cassert>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

double rectangle_trapezoid(function<double(double)> f, double l, double r,
                           double h, bool is_rectangle = true) {
  if (l > r) {
    return 0;
  }
  double result = 0.0;
  double cur_x = l;
  while (cur_x < r) {
    if (is_rectangle) {
      // средних прямоугольников
      result += f((cur_x + cur_x + h) * 0.5);
    } else {
      // трапеций
      result += 0.5 * (f(cur_x + h) + f(cur_x));
    }
    cur_x += h;
  }
  return h * result;
}

double simpson(function<double(double)> f, double l, double r, double h) {
  if (l > r) {
    return 0;
  }
  // корректируем шаг, чтобы интервал делился на четное число частей
  while (static_cast<int>((r - l) / h) % 2 != 0) {
    h *= 0.9;
  }
  double result = 0.0;
  double cur_x = l + h;
  while (cur_x < r) {
    result += f(cur_x - h) + 4.0 * f(cur_x) + f(cur_x + h);
    cur_x += 2 * h;
  }
  return result * h / 3.0;
}

// Рунге-Ромберг-Ричардсон
double RRR(double h1, double h2, double i1, double i2, int p) {
  return i1 + (i1 - i2) / (pow(h2 / h1, p) - 1.0);
}

int main() {
  auto f = [](double x) {
      return x / ((2 * x + 7) * (3 * x + 4));
      //return x/((3*x+4)*(3*x+4));
  };

  ifstream input_file("tests/3-5.json");
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

  double x0 = data["x0"].get<double>();
  double xk = data["xk"].get<double>();
  double h1 = data["h1"].get<double>();
  double h2 = data["h2"].get<double>();
  double p = data["p"].get<double>();

  auto rect1 = rectangle_trapezoid(f, x0, xk, h1);
  auto rect2 = rectangle_trapezoid(f, x0, xk, h2);
  auto trap1 = rectangle_trapezoid(f, x0, xk, h1, false);
  auto trap2 = rectangle_trapezoid(f, x0, xk, h2, false);
  auto s1 = simpson(f, x0, xk, h1);
  auto s2 = simpson(f, x0, xk, h2);

  if (rect1 && rect2 && trap1 && trap2 && s1 && s2) {
    cout << "Rectangle method\n";
    cout << "step " << h1 << ": " << rect1 << endl;
    cout << "step " << h2 << ": " << rect2 << endl;
    cout << "Trapezoid method\n";
    cout << "step " << h1 << ": " << trap1 << endl;
    cout << "step " << h2 << ": " << trap2 << endl;
    cout << "Simpson method\n";
    cout << "step " << h1 << ": " << s1 << endl;
    cout << "step " << h2 << ": " << s2 << endl;

    double r_rrr = RRR(h1, h2, rect1, rect2, p);
    double t_rrr = RRR(h1, h2, trap1, trap2, p);
    double s_rrr = RRR(h1, h2, s1, s2, p);
    cout << "Runge-Rombert-Richardson method\n";
    cout << "Rectangle: " << r_rrr << endl;
    cout << "Trapezoid: " << r_rrr << endl;
    cout << "Simpson: " << r_rrr << endl;
  } else {
    cout << "Integration failed: invalid range" << endl;
  }

  return 0;
}
