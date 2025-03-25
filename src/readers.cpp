#include "../lib/readers.h"

std::pair<Matrix, Matrix> read_SLE(std::ifstream file) {
    if (!file.is_open()) {
        std::cerr << "failed to open a file!" << std::endl;
        exit(0);
    }
    std::size_t rows, cols;
    file >> rows >> cols;
    Matrix A(rows, cols), b(rows, 1);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j <= cols; ++j) {
            if (j < cols) {
                file >> A.data[i][j];
            } else {
                file >> b.data[i][0];
            }
        }
    }
    file.close();
    return std::make_pair(A, b);
}
