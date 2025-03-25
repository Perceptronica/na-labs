#pragma once
#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include "matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(const Matrix &A);
Matrix solveLU(const Matrix &A, const Matrix &b);
Matrix forwardSubstitution(const Matrix &L, const Matrix &b);
Matrix backwardSubstitution(const Matrix &U, const Matrix &y);

Matrix solveTDMA(const Matrix &A, const Matrix &b);

std::pair<Matrix, int> simpleIterationImpl(const Matrix &A, const Matrix &B,
                                           long double epsilon,
                                           int maxIterations, long double tau);
long double vectorDiffNorm(std::vector<long double> &x,
                           std::vector<long double> &y);
std::pair<Matrix, int> Seidel(const Matrix &A, const Matrix &b, long double eps,
                              int maxIters);

#endif
