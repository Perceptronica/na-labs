#include <bits/stdc++.h>

const double EPS = 1e-6;

double f1(double x1, double x2) { return x1 - cos(x2) - 1; }

double f2(double x1, double x2) { return x2 - log10(x1 + 1) - 3; }

void jacobian(double x1, double x2, std::vector<std::vector<double>>& J) {
    J[0][0] = 1.0;                                  // df1/dx1
    J[0][1] = sin(x2);                           // df1/dx2
    J[1][0] = -1.0 / ((x1 + 1) * log(10));       // df2/dx1
    J[1][1] = 1.0;                                  // df2/dx2
}

void newtonMethod(double eps = EPS) {
    double x1 = 1.0, x2 = 3.0;
    std::vector<std::vector<double>> J(2, std::vector<double>(2));
    std::vector<double> F(2);
    std::vector<double> delta(2);
    int iter = 0;
    const int max_iter = 20;

    std::cout << "\nМетод Ньютона:\n";
    std::cout << "iter\tx1\t\tx2\t\tf1\t\tf2\n";

    do {
        jacobian(x1, x2, J);
        F[0] = f1(x1, x2);
        F[1] = f2(x1, x2);

        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (abs(det) < 1e-10) {
            std::cerr << "Якобиан вырожден\n";
            return;
        }

        delta[0] = (F[0] * J[1][1] - F[1] * J[0][1]) / det;
        delta[1] = (J[0][0] * F[1] - J[1][0] * F[0]) / det;

        std::cout << iter << "\t" << std::fixed << std::setprecision(6)
             << x1 << "\t" << x2 << "\t"
             << f1(x1, x2) << "\t" << f2(x1, x2) << std::endl;

        x1 += delta[0];
        x2 += delta[1];

        if (abs(delta[0]) < eps && abs(delta[1]) < eps) {
            break;
        }

        iter++;
    } while (iter < max_iter);

    std::cout << "\nРешение:\nx1 = " << x1 << "\nx2 = " << x2 << std::endl;
    std::cout << "Невязки:\nf1 = " << f1(x1, x2) << "\nf2 = " << f2(x1, x2) << std::endl;
}

void simpleIteration(double eps = EPS) {
    double x1 = 1.0, x2 = 3.0;
    double x1_new, x2_new;
    int iter = 0;
    const int max_iter = 100;

    std::cout << "Метод простой итерации:\n";
    std::cout << "iter\tx1\t\tx2\t\tf1\t\tf2\n";

    do {
        x1_new = cos(x2) + 1;
        x2_new = log10(x1 + 1) + 3;

        std::cout << iter << "\t" << std::fixed << std::setprecision(6)
             << x1 << "\t" << x2 << "\t"
             << f1(x1, x2) << "\t" << f2(x1, x2) << std::endl;

        if (abs(x1_new - x1) < eps && abs(x2_new - x2) < eps) {
            break;
        }

        x1 = x1_new;
        x2 = x2_new;
        iter++;
    } while (iter < max_iter);

    std::cout << "\nРешение:\nx1 = " << x1 << "\nx2 = " << x2 << std::endl;
    std::cout << "Невязки:\nf1 = " << f1(x1, x2) << "\nf2 = " << f2(x1, x2) << std::endl;
}

int main() {
    std::cout << "Решение системы уравнений:\n";
    std::cout << "x1 - cos(x2) = 1\n";
    std::cout << "x2 - lg(x1 + 1) = 3\n\n";

    double eps;
    std::cout << "Введите точность вычислений: ";
    std::cin >> eps;

    simpleIteration(eps);
    newtonMethod(eps);

    return 0;
}
