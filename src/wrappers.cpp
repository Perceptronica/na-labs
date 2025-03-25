#include "../lib/wrappers.h"
#include <stdexcept>

void LU() {
    Matrix A;
    std::cout << "\033[32minput A:\033[0m ";
    std::cin >> A;
    Matrix b(A.rows, 1);
    std::cout << "\033[32minput B:\033[0m ";
    std::cin >> b;
    Matrix x = solveLU(A, b);
    x.transpose();
    std::cout << "\033[33mX = " << x << "^T" << std::endl;
    std::cout << "det(A) = " << A.determinant() << std::endl;
    Matrix C = A.inverse();
    std::cout << "A^(-1) = \n" << C << "\033[0m" << std::endl;
}

void LU(std::pair<Matrix, Matrix> p) {
    Matrix& A = p.first;
    Matrix& b = p.second;
    if (b.cols != 1) {
        throw std::logic_error("b should be a vector");
    }
    std::cout << A << std::endl;
    std::cout << b << std::endl;
    Matrix x = solveLU(A, b);
    x.transpose();
    Matrix C = A.inverse();
    std::cout << x << std::endl;
    std::cout << C << std::endl;
}

void TDMA() {
    Matrix A;
    std::cout << "\033[32minput A:\033[0m ";
    std::cin >> A;
    Matrix b(A.rows, 1);
    std::cout << "\033[32minput B:\033[0m ";
    std::cin >> b;
    Matrix x = solveTDMA(A, b);
    std::cout << "\033[33mX = " << x << "^T" << std::endl;
}

void TDMA(std::pair<Matrix, Matrix> p) {
    Matrix& A = p.first;
    Matrix& b = p.second;
    if (b.cols != 1) {
        throw std::logic_error("b should be a vector");
    }
    std::cout << A << std::endl;
    std::cout << b << std::endl;
    Matrix x = solveTDMA(A, b);
    std::cout << x << std::endl;
}

void SI() {
    Matrix A;
    std::cout << "\033[32minput A:\033[0m ";
    std::cin >> A;
    Matrix b(A.rows, 1);
    std::cout << "\033[32minput B:\033[0m ";
    std::cin >> b;
    std::cout << "enter <eps iterations tau>: ";
    long double e, t;
    int i;
    std::cin >> e >> i >> t;
    auto [x, a] = simpleIterationImpl(A, b, e, i, t);
    std::cout << "iterations: " << a << std::endl;
    std::cout << x << std::endl;
}

void SI(std::pair<Matrix, Matrix> p) {
    Matrix& A = p.first;
    Matrix& b = p.second;
    if (b.cols != 1) {
        throw std::logic_error("b should be a vector");
    }
    std::cout << A << std::endl;
    std::cout << b << std::endl;
    std::cout << "enter <eps iterations tau>: ";
    long double e, t;
    int i;
    std::cin >> e >> i >> t;
    auto [x, a] = simpleIterationImpl(A, b, e, i, t);
    std::cout << "iterations: " << a << std::endl;
    std::cout << x << std::endl;
}

void SeidelW() {
    Matrix A;
    std::cout << "\033[32minput A:\033[0m ";
    std::cin >> A;
    Matrix b(A.rows, 1);
    std::cout << "\033[32minput B:\033[0m ";
    std::cin >> b;
    std::cout << "enter <eps iterations>: ";
    long double e;
    int i;
    std::cin >> e >> i;
    auto [x, a] = Seidel(A, b, e, i);
    std::cout << "iterations: " << a << std::endl;
    std::cout << x << std::endl;
}

void SeidelW(std::pair<Matrix, Matrix> p) {
    Matrix& A = p.first;
    Matrix& b = p.second;
    if (b.cols != 1) {
        throw std::logic_error("b should be a vector");
    }
    std::cout << A << std::endl;
    std::cout << b << std::endl;
    std::cout << "enter <eps iterations>: ";
    long double e;
    int i;
    std::cin >> e >> i;
    auto [x, a] = Seidel(A, b, e, i);
    std::cout << "iterations: " << a << std::endl;
    std::cout << x << std::endl;
}
