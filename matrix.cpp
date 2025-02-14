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
    throw std::invalid_argument("matrix must have at least one element");
  }
}

Matrix::Matrix(std::vector<std::vector<double>> &v) {
  if (v.size() > 0) {
    rows = v.size();
    for (size_t i = 1; i < rows; ++i) {
      if (v[i].size() != v[i - 1].size()) {
        throw std::invalid_argument("matrix must have rows of the same size");
      }
    }
    cols = v[0].size();
    data = v;
  } else {
    throw std::invalid_argument("matrix cannot be empty");
  }
}

Matrix Matrix::operator+(const Matrix &rhs) const {
  if (rows != rhs.rows || cols != rhs.cols) {
    throw std::invalid_argument(
        "matrices must have the same dimensions for addition");
  }
  if (rows == 0 || cols == 0 || rhs.rows == 0 || rhs.cols == 0) {
    throw std::invalid_argument(
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
    throw std::invalid_argument(
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
  if (cols != rhs.rows) {
    throw std::invalid_argument(
        "matrices must have compatible dimensions for multiplication");
  }
  Matrix res(rows, rhs.cols);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < rhs.cols; ++j) {
      double sum = 0;
      for (std::size_t k = 0; k < cols; ++k) {
        sum += data[i][k] * rhs.data[k][j];
      }
      res.data[i][j] = sum;
    }
  }
  return res;
}

Matrix Matrix::operator=(const Matrix &rhs) {
  if (this == &rhs) {
    return *this;
  }
  rows = rhs.rows;
  cols = rhs.cols;
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      data[i][j] = rhs.data[i][j];
    }
  }
  return *this;
}

bool Matrix::operator==(const Matrix &rhs) const {
  if (rows != rhs.rows || cols != rhs.cols) {
    return false;
  }
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      if (data[i][j] != rhs.data[i][j]) {
        return false;
      }
    }
  }
  return true;
}

bool Matrix::operator!=(const Matrix &rhs) const { return !(*this == rhs); }

void Matrix::transpose() {
  Matrix temp(cols, rows);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      temp.data[j][i] = data[i][j];
    }
  }
  *this = temp;
}

double Matrix::determinant() const {
  if (this->rows != this->cols) {
    throw std::invalid_argument("matrix must be square");
  }
  if (this->rows == 1) {
    return data[0][0];
  }
  if (this->rows == 2) {
    return data[0][0] * data[1][1] - data[0][1] * data[1][0];
  }
  double res = 0;
  for (std::size_t i = 0; i < rows; ++i) {
    Matrix submatrix(rows - 1, cols - 1);
    for (std::size_t j = 1; j < rows; ++j) {
      std::size_t k = 0;
      for (std::size_t l = 0; l < cols; ++l) {
        if (l != i) {
          submatrix.data[j - 1][k++] = data[j][l];
        }
      }
    }
    res += data[0][i] * submatrix.determinant() * (i % 2 == 0 ? 1 : -1);
  }
  return res;
}

uint32_t Matrix::rank() const {
  uint32_t rank = 0;
  Matrix temp = *this;
  for (std::size_t i = 0; i < rows; ++i) {
    if (temp.data[i][i] != 0) {
      for (std::size_t j = 0; j < rows; ++j) {
        if (j != i) {
          double factor = temp.data[j][i] / temp.data[i][i];
          for (std::size_t k = i; k < cols; ++k) {
            temp.data[j][k] -= factor * temp.data[i][k];
          }
        }
      }
      ++rank;
    }
  }
  return rank;
}

void Matrix::inverse() {
    // TODO: IMPLEMENT!
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
