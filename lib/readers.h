#include "matrix.h"

#ifndef READ_H
#define READ_H

std::pair<Matrix, Matrix> read_SLE(std::ifstream file);
Matrix read_matrix(std::ifstream file);

#endif
