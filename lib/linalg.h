#pragma once
#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include "matrix.h"

std::pair<Matrix, Matrix> LUDecomposition(const Matrix &A);
Matrix solveLU(const Matrix &A, const Matrix &b);
Matrix forwardSubstitution(const Matrix &L, const Matrix &b);
Matrix backwardSubstitution(const Matrix &U, const Matrix &y);

Matrix solveTDMA(const Matrix &A, const Matrix &b);

long double norm(std::size_t mode, const Matrix &V);
void simpleIterationImpl(const Matrix &A, const Matrix &B, long double epsilon);
void Seidel(const Matrix &a, const Matrix &b, long double epsilon);
void Jacobi(const Matrix &Q, long double epsilon);

#endif
