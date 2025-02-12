#include "matrix.h"

int main() {
    std::vector<std::vector<double>> v = {{1.0, -1.52}, {0.0, 1.0}};
    Matrix a(v);
    Matrix b = a + a + a;
    print(b);
}
