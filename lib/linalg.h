#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include "matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(Matrix& A);
Matrix solveLS(Matrix& A, const std::vector<long double>& b);

#endif
