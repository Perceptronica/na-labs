#include <bits/stdc++.h>
#include <cmath>

double f(double x) {
    return exp(x) - 2*x - 2;
}

double df(double x) {
    return exp(x) - 2;
}

double phi(double x) {
    return log(2*x + 2); //phi'(x) = 1/(1+x) на [1, 2]: |phi'(x)| <= 1/2 = q < 1
}

double newtonMethod(double x0, double eps, int maxIter) {
    double x1;
    int iter = 0;
    do {
        std::cout << std::setw(10) << std::left << iter <<
                "\t" << x0 << "\t" << f(x0) << "\t" << df(x0) << "\t"
                << (-1) * f(x0)/df(x0) << "\n";
        x1 = x0 - f(x0)/df(x0);
        if (fabs(x1 - x0) <= eps) {
            break;
        }
        x0 = x1;
        iter++;
    } while (iter < maxIter);
    if (iter == maxIter) {
        std::cout << "max iters!\n";
    }
    return x1;
}

double simpleIters(double k, double l, double q, double eps, int maxIter) {
    double x1;
    int iter = 0;
    double x0 = (k + l) / 2.0;
    do {
        std::cout << std::setw(10) << std::left << iter <<
                "\t" << x0 << "\t" << phi(x0) << "\n";
        x1 = phi(x0);
        if ((q/(1-q)) * fabs(x1 - x0) <= eps) {
            break;
        }
        x0 = x1;
        iter++;
    } while (iter < maxIter);
    if (iter == maxIter) {
        std::cout << "max iters!\n";
    }
    return x1;
}

int main() {
    double x0 = 2.0;    //f(2.0)*f''(2.0) > 0
    double eps = 1e-4;
    int maxIter = 20;
    double a = 1;
    double b = 2;
    double q = 0.5;
    std::cout << "Newthon method:\n";
    double ans1 = newtonMethod(x0, eps, maxIter);
    std::cout << "Simple iterations method:\n";
    double ans2 = simpleIters(a, b, q, eps, maxIter);
    return 0;
}
