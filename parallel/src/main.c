//
// Created by evgenii on 27.11.17.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
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

void generate_matrix_vector(int n, double** A, double** B) {
    double* Start_A = calloc((size_t) ((unsigned long) n * (unsigned long) n), sizeof(double));
    double* Start_B = calloc((size_t) n, sizeof(double));
    for (int i = 0; i < n; ++i) {
        Start_B[i] = getrand(2, 700);
        for (int j = 0; j < n; ++j) {
            if (i == j)
                Start_A[i * n + j] = 100;
            else {
                int r = 0;
                while ((r = getrand(-40, 40)) == 0);
                Start_A[i * n + j] = r;
            }
        }
    }
    *A = Start_A;
    *B = Start_B;
}

int main(int argc, char** argv) {
    double time_total = -wtime();
    MPI_Init(&argc, &argv);
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int N ,M;
    double eps;
    if (rank == 0) {
        if (argc < 4) {
            perror("Error: Wrong number of arguments.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        N = (int) strtol(argv[1], NULL, 10);
        M = (int) strtol(argv[2], NULL, 10);
        if (N != M) {
            perror("Error: N != M");
        }
        eps = pow(10, -strtol(argv[3], NULL, 10));
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int* nrows = calloc((size_t) commsize, sizeof(int));
    int* displs = calloc((size_t) commsize, sizeof(int));
    for (int i = 0; i < commsize; ++i) {
        int l, u;
        get_chunk(0, N - 1, commsize, i, &l, &u);
        nrows[i] = u - l + 1;
        displs[i] = (i > 0) ? displs[i - 1] + nrows[i - 1] : 0;
    }

    double time_init = -wtime();
    double *A = calloc((size_t) nrows[rank] * M, sizeof(double));
    double *B = calloc((size_t) nrows[rank], sizeof(double));
    double *X = calloc((size_t) M, sizeof(double));
    time_init += wtime();

    double time_fill = -wtime();
    for (int i = 0; i < M; ++i) {
        X[i] = 7;
    }
    if (rank == 0) {
        FILE* fd_m;
        FILE* fd_v;
        unsigned int global_iter_A = 0;
        unsigned int global_iter_B = 0;
        double* Start_A;
        double* Start_B;
        if (argc == 6) {
            fd_m = fopen(argv[4], "r");
            fd_v = fopen(argv[5], "r");
            if (!fd_m || !fd_v) {
                perror("Can't open file!");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        else if (argc == 4) {
            generate_matrix_vector(N, &Start_A, &Start_B);
        }
        else {
            perror("Wrong number of arguments!");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < nrows[0]; ++i) {
            if (argc == 6)
                fscanf(fd_v, "%lf", &B[i]);
            else if (argc == 4) {
                B[i] = Start_B[global_iter_B++];
            }
            for (int j = 0; j < M; ++j) {
                if (argc == 6)
                    fscanf(fd_m, "%lf", &A[i * M + j]);
                else if (argc == 4) {
                    A[i * M + j] = Start_A[global_iter_A++];
                }
            }
        }
        for (int r = 1; r < commsize; ++r) {
            int size = nrows[r] * M;
            double *PART_A = calloc((size_t) size, sizeof(double));
            double *PART_B = calloc((size_t) nrows[r], sizeof(double));
            for (int i = 0; i < nrows[r]; ++i) {
                if (argc == 6)
                    fscanf(fd_v, "%lf", &PART_B[i]);
                else if (argc == 4) {
                    PART_B[i] = Start_B[global_iter_B++];
                }
                for (int j = 0; j < M; ++j) {
                    if (argc == 6)
                        fscanf(fd_m, "%lf", &PART_A[i * M + j]);
                    else if (argc == 4) {
                        PART_A[i * M + j] = Start_A[global_iter_A++];
                    }
                }
            }
            MPI_Send(PART_A, size, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
            MPI_Send(PART_B, nrows[r], MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
            free(PART_A);
            free(PART_B);
        }
        if (argc == 6) {
            fclose(fd_m);
            fclose(fd_v);
        }
        else if (argc == 4) {
            free(Start_A);
            free(Start_B);
        }
    }
    else {
        MPI_Recv(A, (nrows[rank] * M), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(B, nrows[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    time_fill += wtime();
    double time_jacobi = -wtime();
    Jacobi(M, A, B, X, eps, nrows, displs);
    time_jacobi += wtime();
    time_total += wtime();
    if (rank == 0) {
        printf("Jacobi methood (%dX%d) proc %d\n"
                       "Time Total: %.6f\n"
                       "Time Jacobi %.6f\n"
                       "Time Init %.6f\n"
                       "Time Fill %.6f\n", N, M, commsize, time_total, time_jacobi, time_init, time_fill);
    }
    free(displs);
    free(nrows);
    free(A);
    free(B);
    free(X);
    MPI_Finalize();
    return 0;
}