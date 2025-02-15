#pragma once
#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include "matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(const Matrix& A);
Matrix solveLU(const Matrix& A, const Matrix& b);

#endif
