#include "lib/wrappers.h"
#include "lib/readers.h"

int main(int argc, char *argv[]) {
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
      if (file_path.empty()) { LU(); }
      else {
          LU(read_SLE(std::ifstream(file_path)));
      }
    } else if (command == "exit") {
      exit = true;
    }
  }
}
