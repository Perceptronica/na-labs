#include "lib/readers.h"
#include "lib/wrappers.h"

int main() {
  std::cout << "\033[32mNumerical Analysis labs\nby Anton Kabanov "
               "(@Perceptronica)\033[0m"
            << std::endl;
  std::string command;
  std::string file_path;
  std::string input;
  bool exit = false;
  while (!exit) {
    std::cout << "> ";
    std::getline(std::cin, input);
    std::istringstream iss(input);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }
    if (tokens.empty()) {
      continue;
    }
    std::string command = tokens[0];
    std::string file_path;
    if (tokens.size() > 1) {
      file_path = tokens[1];
    }

    if (command == "lu") {
      if (file_path.empty()) {
        LU();
      } else {
        LU(read_SLE(std::ifstream(file_path)));
      }
    } else if (command == "td" || command == "tdma") {
      if (file_path.empty()) {
        TDMA();
      } else {
        TDMA(read_SLE(std::ifstream(file_path)));
      }
    } else if (command == "si") {
      if (file_path.empty()) {
        SI();
      } else {
        SI(read_SLE(std::ifstream(file_path)));
      }
    } else if (command == "seidel") {
      if (file_path.empty()) {
        SeidelW();
      } else {
        SeidelW(read_SLE(std::ifstream(file_path)));
      }
    } else if (command == "jacobi") {
      if (file_path.empty()) {
        Matrix S;
        std::cin >> S;
        long double e;
        std::cout << "enter epsilon: ";
        std::cin >> e;
        Jacobi(S, e);
      } else {
        Matrix M = read_matrix(std::ifstream(file_path));
        //std::cout << M;
        long double e;
        std::cout << "enter epsilon: ";
        std::cin >> e;
        Jacobi(M, e);
      }
    } else if (command == "exit") {
      exit = true;
    }
  }
}
