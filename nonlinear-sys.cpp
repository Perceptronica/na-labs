#include "norms.h"
#include <bits/stdc++.h>

using namespace std;

const double EPS = 1e-6;

double f1(double x1, double x2) { return x1 - cos(x2) - 1; }

double f2(double x1, double x2) { return x2 - log10(x1 + 1) - 3; }

void jacobian(double x1, double x2, std::vector<std::vector<long double>> &J) {
  J[0][0] = 1.0;                         // df1/dx1
  J[0][1] = sin(x2);                     // df1/dx2
  J[1][0] = -1.0 / ((x1 + 1) * log(10)); // df2/dx1
  J[1][1] = 1.0;                         // df2/dx2
}

void simpleIteration(double eps = EPS) {
  double x1 = 0, x2 = 3.5;
  double x1_new, x2_new;
  int iter = 0;
  const int max_iter = 100;

  std::cout << "Метод простой итерации:\n";
  std::cout << "iter\tx1\t\tx2\t\tf1\t\tf2\n";

  do {
    x1_new = cos(x2) + 1;
    x2_new = log10(x1 + 1) + 3;

    std::cout << iter << "\t" << std::fixed << std::setprecision(6) << x1
              << "\t" << x2 << "\t" << f1(x1, x2) << "\t" << f2(x1, x2)
              << std::endl;

    if (abs(x1_new - x1) <= eps && abs(x2_new - x2) <= eps) {
      break;
    }

    x1 = x1_new;
    x2 = x2_new;
    iter++;
  } while (iter < max_iter);

  std::cout << "\nРешение:\nx1 = " << x1 << "\nx2 = " << x2 << std::endl;
  std::cout << "f(X):\nf1 = " << f1(x1, x2) << "\nf2 = " << f2(x1, x2)
            << std::endl;
}

void Newton() {
  auto f1 = [](double x1, double x2) { return x1 - cos(x2) - 1; };
  auto f2 = [](double x1, double x2) { return x2 - log10(x1 + 1) - 3; };
  double x1 = 1.0;
  double x2 = 1.0;
  const double tol = 1e-6;
  const int max_iter = 50;
  const double ln10 = log(10.0);

  cout << "iter\t   x1\t\t   x2\t\t   f1\t\t   f2" << endl;
  cout << string(60, '-') << endl;
  cout << fixed << setprecision(7);

  for (int i = 0; i < max_iter; ++i) {
    double func1 = f1(x1, x2);
    double func2 = f2(x1, x2);
    double norm = sqrt(func1 * func1 + func2 * func2);

    cout << setw(2) << i << "\t" << setw(12) << x1 << "\t" << setw(12) << x2
         << "\t" << setw(12) << func1 << "\t" << setw(12) << func2 << endl;

    if (norm < tol) {
      cout << "x1 = " << x1 << ", x2 = " << x2 << endl;
      return;
    }
    double j11 = 1.0;                      // df1/dx1
    double j12 = sin(x2);                  // df1/dx2
    double j21 = -1.0 / ((x1 + 1) * ln10); // df2/dx1
    double j22 = 1.0;                      // df2/dx2
    double det = j11 * j22 - j12 * j21;

    double dx1 = (-func1 * j22 - (-func2) * j12) / det;
    double dx2 = (j11 * (-func2) - j21 * (-func1)) / det;
    x1 += dx1;
    x2 += dx2;
  }
  cout << "\nМаксимальное число итераций достигнуто" << endl;
}

int main() {
  std::cout << "Решение системы уравнений:\n";
  std::cout << "x1 - cos(x2) = 1\n";
  std::cout << "x2 - lg(x1 + 1) = 3\n";

  double eps;
  std::cout << "Введите точность вычислений: ";
  std::cin >> eps;
  // (0.00944, 3.00408)
  simpleIteration(eps);
  Newton();
  return 0;
}
