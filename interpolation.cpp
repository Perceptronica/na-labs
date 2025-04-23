#include <bits/stdc++.h>
#include <cmath>

long double f(long double x) { return exp(x); }

std::vector<long double> X_1 = {-2, -1, 0, 1};
std::vector<long double> X_2 = {-2, -1, 0.2, 1};
long double X_star = -0.5;

long double Lagrangia(std::vector<std::vector<long double>>& table, long double X_star){
    std::size_t dim = table[0].size();
    std::size_t n = dim - 1;
    table.push_back({});
    for (std::size_t i = 0; i <= n; ++i) { // w'
        long double res = 1.0;
        for (std::size_t j = 0; j <= n; ++j) {
            if (i != j) {
                res *= (table[0][i] - table[0][j]);
            }
        }
        table[2].push_back(res);
    }
    table.push_back({});
    for (std::size_t i = 0; i <= n; ++i) { // f / w'
        long double res = table[1][i] / table[2][i];
        table[3].push_back(res);
    }
    table.push_back({});
    for (std::size_t i = 0; i <= n; ++i) { // X* - x_i
        long double res = X_star - table[0][i];
        table[4].push_back(res);
    }
    for (std::size_t i = 0; i < table.size() - 1; ++i) {
        for (std::size_t j = 0; j <= dim; ++j) {
            std::cout << std::setw(12) << std::left << table[j][i];
        }
        std::cout << "\n";
    }
    std::cout << "L_" << n << "(x) = ";
    for (std::size_t i = 0; i <= n; ++i) {
        std::cout << table[3][i] << "*";
        for (std::size_t j = 0; j <= n; ++j) {
            if (i != j) {
                std::cout << "(x+" << table[0][j] << ")";
            }
        }
        if (i != n) {
            std::cout << " + ";
        }
    }
    std::cout << "\n";
    long double res = 0.0;
    for (std::size_t i = 0; i <= n; ++i) {
        long double q = table[3][i];
        long double w = 1.0;
        for (std::size_t j = 0; j <= n; ++j) {
            if (i != j) {
                w *= (X_star - table[0][j]);
            }
        }
        res += q * w;
    }
    std::cout << "L_" << n << "(" << X_star << ") = " << res << "\n";
    return res;
}

long double Newton(std::vector<std::vector<long double>>& table, long double X_star) {
    //...
}

int main() {
    std::vector<std::vector<long double>> t;
    t.push_back(X_1);
    t.push_back({});
    for (std::size_t i = 0; i < X_1.size(); ++i) {
        t[1].push_back(f(X_1[i]));
    }
    double x = Lagrangia(t, X_star);
    double y = f(X_star);
    std::cout << "y(" << X_star << ") = " << y << "\n";
    std::cout << "Абсолютная погрешность интерполяции: " << fabs(y - x) << "\n";
}
