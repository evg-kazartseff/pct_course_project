//
// Created by evgenii on 27.11.17.
//

#ifndef MPI_SIMPLE_ITERATION_JACOBI_H
#define MPI_SIMPLE_ITERATION_JACOBI_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <zconf.h>
#include <mpi.h>

#define PRINT 0

void get_chunk(int a, int b, int commsize, int rank, int* lb, int* ub);
void Jacobi(int M, const double* A, const double* B, double* X, double eps, int* nrows, int* displs);


#endif //MPI_SIMPLE_ITERATION_JACOBI_H
