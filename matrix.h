#pragma once
#include <bits/stdc++.h>
#include <iterator>
#include <stdexcept>

#ifndef _MATRIX_H_
#define _MATRIX_H_

struct Matrix {
  std::size_t rows;
  std::size_t cols;
  std::vector<std::vector<long double>> data;

  Matrix() : rows(0), cols(0), data(0) {};
  Matrix(std::size_t r, std::size_t c);
  Matrix(std::vector<long double> &v); // for vectors/rows
  Matrix(std::vector<std::vector<long double>> &v);
  Matrix(const Matrix &other);
  Matrix(Matrix &&other) noexcept;

  ~Matrix() = default;

  Matrix operator+(const Matrix &rhs) const;
  Matrix operator-(const Matrix &rhs) const;
  Matrix operator*(const double &t) const;
  Matrix operator*(const Matrix &rhs) const;
  Matrix &operator=(const Matrix &other);
  bool operator==(const Matrix &rhs) const;
  bool operator!=(const Matrix &rhs) const;

  friend std::istream &operator>>(std::istream &is, Matrix &mat);
  friend std::ostream &operator<<(std::ostream &os, const Matrix &mat);

  double determinant() const;
  void transpose();
  uint32_t rank() const;
  Matrix inverse();

  Matrix getMinor(std::size_t row, std::size_t col) const;
  Matrix cofactorMatrix() const;
  Matrix adjugateMatrix() const;
};

#endif
