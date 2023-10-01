/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* durbin.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/* Include polybench common header. */
#include "polybench.h"

/* Include benchmark-specific header. */
#include "durbin.h"


/* Variáveis que tínhamos no arquivo header antes de passar o tamanho do dataset como argumento de CLI */
int N;
int _PB_N;
int world_size;
int world_rank;


/* Mensagem de ajuda do programa (--help) */
void help_message() {
    printf("--- Instruções de uso ---\n");
    printf("./durbin.out [OPÇÕES]\n");
    printf("\n[OPÇÕES]\n");
    printf("-d: Tamanho do dataset. Pode ser um dentre {small, medium, large}.\n");
}

/* Array initialization. */
static void init_array(int n, DATA_TYPE POLYBENCH_1D(r,N,n)) {
    int i, j;

    for (i = 0; i < n; i++) {
        r[i] = (n+1-i);
    }
}

/* DCE code. Must scan the entire live-out data.
 *Can be used also to check the correctness of the output. */
static void print_array(int n, DATA_TYPE POLYBENCH_1D(y,N,n)) {
    int i;

    for (i = 0; i < n; i++) {
        printf(DATA_PRINTF_MODIFIER, y[i]);
    }

    printf("\n");
}


/* Main computational kernel. The whole function will be timed,
 *including the call and return. */
static void kernel_durbin(int n,
                          DATA_TYPE POLYBENCH_1D(r,N,n),
                          DATA_TYPE POLYBENCH_1D(y,N,n)) {
    DATA_TYPE z[N];
    DATA_TYPE alpha;
    DATA_TYPE beta;
    DATA_TYPE sum;
    int i, k;
    DATA_TYPE buf;

    #pragma scop
    y[0] = -r[0];
    beta = SCALAR_VAL(1.0);
    alpha = -r[0];

    for (k = 1; k < _PB_N; k++) {
        beta = (1-alpha*alpha)*beta;
        sum = SCALAR_VAL(0.0);

        MPI_Barrier(MPI_COMM_WORLD);

        for (i = world_rank * k / world_size; i < (k + k * world_rank) / world_size; i++) {
            sum += r[k-i-1]*y[i];  // condição de corrida - lock
        }

        // cada processo envia sua soma parcial para o processo 0 fazer o somatório e fazer o broadcast do valor final para todos
        if (world_rank == 0) {
            for (int j = 0; j < world_size - 1; j++) {
                MPI_Recv(&buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sum += buf;
            }
        } else {
            MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        alpha = - (r[k] + sum)/beta;

        for (i = world_rank * k / world_size; i < (k + k * world_rank) / world_size; i++) {
            z[i] = y[i] + alpha*y[k-i-1];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        for (i = world_rank * k / world_size; i < (k + k * world_rank) / world_size; i++) {
            y[i] = z[i];
        }

        y[k] = alpha;

        // if p0: recebo a parte de cada processo do vetor y e escrevo nas posições corretas; envio pra todo mundo vetor y completo
        // if nao p0: send pra p0 da minha parte de y; recebo vetor y inteiro de p0
        if (world_rank == 0) {
            for (int j = 1; j < world_size; j++) {
                for (i = j * k / world_size; i < (k + k * j) / world_size; i++) {
                    MPI_Recv(&y[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // to-do: ver como que faz pra mandar o vetor todo
        } else {
            for (i = world_rank * k / world_size; i < (k + k * world_rank) / world_size; i++) {
                MPI_Send(&y[i], 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
            }

            MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    #pragma endscop
}


int main(int argc, char** argv) {
    /*Lidar com argumentos passados por CLI*/
    /*Lidar com argumentos passados por CLI*/
    for (int i = 0; i < argc; i++) {
        /* Define tamanho da entrada */
        if (strcmp(argv[i], "-d") == 0) {
            i++;

            if (strcmp(argv[i], "small") == 0) {
                N = 16;  // 250k antes
            } else if (strcmp(argv[i], "medium") == 0) {
                N = 350000;
            } else if (strcmp(argv[i], "large") == 0) {
                N = 450000;
            } else {
                help_message();

                return -1;
            }
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) { /* Informações help */
            help_message();

            return 0;
        }
    }

    /* Retrieve problem size. */
    int n = N;
    _PB_N = POLYBENCH_LOOP_BOUND(N,n);

    // inicializações MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    /* Variable declaration/allocation. */
    POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);
    POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);

    /* Initialize array(s). */
    init_array (n, POLYBENCH_ARRAY(r));

    /* Start timer. */
    polybench_start_instruments;

    /* Run kernel. */
    kernel_durbin (n,
                   POLYBENCH_ARRAY(r),
                   POLYBENCH_ARRAY(y));

    /* Stop and print timer. */
    polybench_stop_instruments;
    polybench_print_instruments;

    if (world_rank == 0) {
        print_array(n, POLYBENCH_ARRAY(y));
    }

    /* Be clean. */
    POLYBENCH_FREE_ARRAY(r);
    POLYBENCH_FREE_ARRAY(y);

    MPI_Finalize();

    return 0;
}
