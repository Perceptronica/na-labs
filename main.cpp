#include "lib/wrappers.h"

int main(int argc, char *argv[]) {
  std::cout << "\033[32mNumerical Analysis labs\nby Anton Kabanov "
               "(@Perceptronica)\033[0m"
            << std::endl;
  std::string command;
  bool exit = false;
  while (!exit) {
    std::cout << "> ";
    std::cin >> command;
    if (command == "lu") {
      LU();
    } else if (command == "exit") {
      exit = true;
    }
  }
}
