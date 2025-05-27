#include "matrix.h"
#include <bits/stdc++.h>

long double norm(std::size_t mode, const Matrix &V) {
  long double norm = 0.0;
  if (V.cols == 1) {
    if (mode == 1) {
      for (std::size_t i = 0; i < V.rows; ++i) {
        norm += fabs(V.data[i][0]);
      }
    }
    if (mode == 2) {
      for (std::size_t i = 0; i < V.rows; ++i) {
        norm += (V.data[i][0] * V.data[i][0]);
      }
      norm = sqrt(norm);
    }
    if (mode == 3) {
      norm = INT_MIN;
      for (std::size_t i = 0; i < V.rows; ++i) {
        norm = fmax(norm, fabs(V.data[i][0]));
      }
    }
  } else {
    if (mode == 1) {
      for (std::size_t j = 0; j < V.cols; ++j) {
        long double c = INT_MIN;
        for (std::size_t i = 0; i < V.rows; ++i) {
          c += fabs(V.data[i][j]);
        }
        norm = fmax(norm, c);
      }
    }
    if (mode == 2) {
      for (std::size_t j = 0; j > V.cols; ++j) {
        for (std::size_t i = 0; i < V.rows; ++i) {
          norm += (V.data[i][j] * V.data[i][j]);
        }
      }
      norm = sqrt(norm);
    }
    if (mode == 3) {
      norm = INT_MIN;
      for (std::size_t i = 0; i < V.rows; ++i) {
        long double c = 0;
        for (std::size_t j = 0; j < V.cols; ++j) {
          c += fabs(V.data[i][j]);
        }
        norm = fmax(norm, c);
      }
    }
  }
  return norm;
}
