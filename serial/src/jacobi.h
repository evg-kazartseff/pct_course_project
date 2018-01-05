//
// Created by evgenii on 27.11.17.
//

#ifndef MPI_SIMPLE_ITERATION_JACOBI_H
#define MPI_SIMPLE_ITERATION_JACOBI_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define PRINT_ITER 0

void Jacobi(long N, long M, const double** A, const double* B, double* X, double eps);


#endif //MPI_SIMPLE_ITERATION_JACOBI_H
