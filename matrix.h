#pragma once
#include <bits/stdc++.h>
#include <stdexcept>

#ifndef _MATRIX_H_
#define _MATRIX_H_

struct Matrix {
  std::size_t rows;
  std::size_t cols;
  std::vector<std::vector<double>> data;

  Matrix();
  Matrix(std::size_t r = 1, std::size_t c = 1);
  Matrix(std::vector<std::vector<double>>& v);

  Matrix operator+(const Matrix& rhs) const;
  Matrix operator-(const Matrix& rhs) const;
  Matrix operator*(const double& t) const;
  Matrix operator*(const Matrix& rhs) const;
  Matrix operator=(const Matrix& rhs);
  bool operator==(const Matrix& rhs) const;
  bool operator!=(const Matrix& rhs) const;

  double determinant() const;
  void transpose();
  uint32_t rank() const;
  void inverse();
};

void print(const Matrix& m);

#endif
