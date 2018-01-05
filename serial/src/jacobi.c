//
// Created by evgenii on 27.11.17.
//
#include <zconf.h>
#include "jacobi.h"

void Jacobi(long N, long M, const double** A, const double* B, double* X, double eps)
{
    double* TempX = malloc(N * sizeof(double));
    double norm;
    char format[80];
    sprintf(format, "%%.%.flf\n", -log10(eps));
    sync();
    do {
        for (int i = 0; i < N; i++) {
            TempX[i] = B[i];
            for (int g = 0; g < M; g++) {
                if (i != g)
                    TempX[i] -= A[i][g] * X[g];
            }
            TempX[i] /= A[i][i];
        }
        norm = fabs(X[0] - TempX[0]);
        for (int h = 0; h < N; h++) { //норму не могут считать все процессы, толькоо один
            if (fabs(X[h] - TempX[h]) > norm)
                norm = fabs(X[h] - TempX[h]);
            X[h] = TempX[h];
#if PRINT_ITER
            printf(format, X[h]);
#endif
        }
#if PRINT_ITER
        printf("\n");
        sync();
#endif
    } while (norm  >= eps);
    for (int i = 0; i < N; ++i) {
        printf(format, X[i]);
    }
    printf("\n");
    free(TempX);
}
