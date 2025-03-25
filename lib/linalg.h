#pragma once
#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include "matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(const Matrix& A);
Matrix solveLU(const Matrix& A, const Matrix& b);
Matrix forwardSubstitution(const Matrix &L, const Matrix &b);
Matrix backwardSubstitution(const Matrix &U, const Matrix &y);

#endif
