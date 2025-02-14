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
      std::cout << "A: ";
      std::cin >> A;
      std::vector<long double> b(A.rows);
      std::cout << "b: ";
      for (std::size_t i = 0; i < A.rows; ++i) {
        std::cin >> b[i];
      }
      Matrix x = solveLS(A, b);
      x.transpose();
      std::cout << "\033[34mx:\n" << x << "\033[0m";
      std::cout << "det(A) = " << A.determinant() << std::endl;
    } else if (command == "exit") {
      exit = true;
    }
  }
}
