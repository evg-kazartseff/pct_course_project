//
// Created by evgenii on 27.11.17.
//
#include "jacobi.h"

void get_chunk(int a, int b, int commsize, int rank, int* lb, int* ub)
{
    int n = b - a + 1;
    int q = n / commsize;
    if (n % commsize)
        q++;
    int r = commsize * q - n;
    /* Compute chunk size for the process */
    int chunk = q;
    if (rank >= commsize - r) chunk = q - 1;
    *lb = a; /* Determine start item for the process */
    if (rank > 0) { /* Count sum of previous chunks */
        if (rank <= commsize - r)
            *lb += q * rank;
        else
            *lb += q * (commsize - r) + (q - 1) * (rank - (commsize - r));
    }
    *ub = *lb + chunk - 1;
}

void Jacobi(int M, const double* A, const double* B, double* X, double eps, int* nrows, int* displs)
{
    int rank, commsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    int N = nrows[rank];
    double* TempX = malloc(N * sizeof(double));
    char format[80];
    sprintf(format, "%%.%.flf\n", -log10(eps));
    double norm, global_norm;
    do {
        for (int i = 0; i < N; i++) {
            TempX[i] = B[i];
            for (int g = 0; g < M; g++) {
                if (displs[rank] + i != g) {
                    TempX[i] -= A[i * M + g] * X[g];
                }
            }
            TempX[i] /= A[(i * M + i) + displs[rank]];
        }
        norm = fabs(X[displs[rank]] - TempX[0]);
        for (int i = 0; i < N; i++) {
            if (fabs(X[displs[rank] + i] - TempX[i]) > norm)
                norm = fabs(X[displs[rank] + i] - TempX[i]);
        }
        MPI_Allgatherv(TempX, N, MPI_DOUBLE, X, nrows, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allreduce((float* ) &norm, (float* ) &global_norm, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    } while (global_norm > eps);
#if PRINT
    if (rank == 0) {
        for (int i = 0; i < M; ++i) {
            printf(format, X[i]);
        }
        printf("\n");
    }
#endif
    free(TempX);
}
