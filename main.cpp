#include "lib/linalg.h"
#include "lib/matrix.h"

int main() {
  std::cout << "\033[32mNumerical Analysis labs\nby Anton Kabanov (@Perceptronica)\033[0m" << std::endl;
  std::string command;
  bool exit = false;
  while (!exit) {
    std::cout << "> ";
    std::cin >> command;
    if (command == "lu") {
      std::cout << "LU decomposition | Solution to Ax = b" << std::endl;
      Matrix A;
      std::cout << "\033[32minput A:\033[0m ";
      std::cin >> A;
      Matrix b(A.rows, 1);
      std::cout << "\033[32minput B:\033[0m ";
      std::cin >> b;
      Matrix x = solveLU(A, b);
      x.transpose();
      std::cout << "\033[32mX = " << x << "^T" << std::endl;
      std::cout << "det(A) = " << A.determinant() << "\033[0m"<< std::endl;
      //A.inverse();
      //std::cout << "A^(-1) = \n" << A << std::endl;
    } else if (command == "exit") {
      exit = true;
    }
  }
}
