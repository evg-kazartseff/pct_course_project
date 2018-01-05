//
// Created by evgenii on 27.11.17.
//
//#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "jacobi.h"

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

int getrand(int min, int max)
{
    return (int) (rand() / (RAND_MAX + 1.0) * (max - min) + min);
}

int Fill_matrix(double** Matrix, long N, long M, const char* file) {
    FILE* fd = fopen(file, "r");
    if (!fd) {
        perror("Can't open file!");
        return 1;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            fscanf(fd, "%lf", &Matrix[i][j]);
        }
    }
    fclose(fd);
    return 0;
}

int Fill_B(double* B, long N, const char* file) {
    FILE* fd = fopen(file, "r");
    if (!fd) {
        perror("Can't open file!");
        return 1;
    }
    for (int i = 0; i < N; ++i) {
            fscanf(fd, "%lf", &B[i]);
    }
    fclose(fd);
    return 0;
}

void generate_matrix_vector() {
    FILE* fdm = fopen("10X10/matrix.txt", "w");
    FILE* fdv = fopen("10X10/B.txt", "w");
    int n = 10;
    for (int i = 0; i < n; ++i) {
        fprintf(fdv, "%d\n", getrand(2, 30));
        for (int j = 0; j < n; ++j) {
            if (i == j)
                fprintf(fdm, "%d ", 100);
            else {
                int i = 0;
                while ((i = getrand(-6, 6)) == 0);
                    fprintf(fdm, "%d ", i);

            }

        }
        fprintf(fdm, "\n");
    }
    //fflush(fdm);
    //fflush(fdv);
    fclose(fdm);
    fclose(fdv);
}

int main(int argc, char** argv) {
    srand((unsigned int) time(NULL));
    //generate_matrix_vector();
#if 1
    if (argc != 6) {
        perror("Error: Wrong number of arguments.\n");
        return 1;
    }
    long N = strtol(argv[1], NULL, 10);
    long M = strtol(argv[2], NULL, 10);
    double eps = pow(10, -strtol(argv[3], NULL, 10));
    double time_total = -wtime();
    double time_init = -wtime();
    const double** A = malloc(N * sizeof(double*));
    double *B = calloc((size_t) M, sizeof(double));
    double *X = calloc((size_t) N, sizeof(double));
    for (int i = 0; i < N; ++i) {
        A[i] = calloc((size_t) M, sizeof(double));
    }
    time_init += wtime();
    double time_fill = -wtime();
    for (int i = 0; i < N; ++i) {
        X[i] = 7;
    }
    Fill_matrix(A, N, M, argv[4]);
    Fill_B(B, M, argv[5]);
    time_fill += wtime();
    double time_jacobi = -wtime();
    Jacobi(N, M, A, B, X, eps);
    time_jacobi += wtime();
    time_total +=wtime();
    printf( "Time Total: %.5f\n"
            "Time Jacobi %.5f\n"
            "Time Init %.5f\n"
            "Time Fill %.5f\n", time_total, time_jacobi, time_init, time_fill);
#endif
    return 0;
}