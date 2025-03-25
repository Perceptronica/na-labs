#include "linalg.h"
#include "matrix.h"

#ifndef WRAPPERS_H
#define WRAPPERS_H

void LU();
void LU(std::pair<Matrix, Matrix> p);

void TDMA();
void TDMA(std::pair<Matrix, Matrix> p);

#endif
