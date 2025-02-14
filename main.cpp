#include "matrix.h"

int main() {
    std::vector<std::vector<double>> v = {{1, 2, 3}, {4, 3, 6}, {7, 8, 9}};
    Matrix a(v);
    double det = a.determinant();
    uint32_t rank = a.rank();
    a.transpose();
    std::cout << det << ' ' << rank << std::endl;
    print(a);
}
