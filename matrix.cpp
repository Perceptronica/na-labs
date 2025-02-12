#include "matrix.h"
#include <stdexcept>

Matrix::Matrix() {
  rows = 0;
  cols = 0;
}

Matrix::Matrix(std::size_t r, std::size_t c) {
  if (r > 0 && c > 0) {
    rows = r;
    cols = c;
    data.resize(r, std::vector<double>(c, 0.0));
  } else {
    throw std::logic_error("matrix must have at least one element");
  }
}

Matrix::Matrix(std::vector<std::vector<double>> &v) {
  if (v.size() > 0) {
    rows = v.size();
    for (size_t i = 1; i < rows; ++i) {
      if (v[i].size() != v[i - 1].size()) {
        throw std::logic_error("matrix must have rows of the same size");
      }
    }
    cols = v[0].size();
    data = v;
  } else {
    throw std::logic_error("matrix cannot be empty");
  }
}

Matrix Matrix::operator+(const Matrix &rhs) const {
  if (rows != rhs.rows || cols != rhs.cols) {
    throw std::logic_error(
        "matrices must have the same dimensions for addition");
  }
  if (rows == 0 || cols == 0 || rhs.rows == 0 || rhs.cols == 0) {
    throw std::logic_error(
        "the addition cannot be performed on the empty matrices");
  }
  Matrix res(rows, cols);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      res.data[i][j] = data[i][j] + rhs.data[i][j];
    }
  }
  return res;
}

Matrix Matrix::operator*(const double &t) const {
  if (rows == 0 || cols == 0) {
    throw std::logic_error(
        "the multiplication cannot be performed on an empty matrix");
  }
  Matrix res(rows, cols);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      res.data[i][j] = data[i][j] * t;
    }
  }
  return res;
}

Matrix Matrix::operator-(const Matrix &rhs) const {
  return *this + (rhs * (-1));
}

Matrix Matrix::operator*(const Matrix &rhs) const {
  // ! IMPLEMENT
  return *this + rhs;
}

void print(const Matrix &m) {
  if (m.rows == 0 || m.cols == 0) {
    std::cout << "empty matrix" << std::endl;
  } else {
    for (std::size_t i = 0; i < m.rows; ++i) {
      for (std::size_t j = 0; j < m.cols; ++j) {
        std::cout << m.data[i][j] << '\t';
      }
      std::cout << '\n';
    }
  }
}
