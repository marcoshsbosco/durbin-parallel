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

/* Include polybench common header. */
#include "polybench.h"

/* Include benchmark-specific header. */
#include "durbin.h"


/* Variáveis que tínhamos no arquivo header antes de passar o tamanho do dataset como argumento de CLI */
int N;
int _PB_N;


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

    #pragma scop
    y[0] = -r[0];
    beta = SCALAR_VAL(1.0);
    alpha = -r[0];

    // prints para debug
    printf("\n-- Antes de tudo --\n");
    printf("alpha: %lf\n", alpha);
    printf("beta: %lf\n", beta);
    printf("sum: %lf\n", sum);

    printf("r: ");
    for (int i = 0; i < N; i++) {
        printf("%0.2lf ", r[i]);
    }
    printf("\n");

    printf("y: ");
    for (int i = 0; i < N; i++) {
        printf("%0.2lf ", y[i]);
    }
    printf("\n");
    // fim prints para debug

    for (k = 1; k < _PB_N; k++) {
        beta = (1-alpha*alpha)*beta;
        sum = SCALAR_VAL(0.0);

        printf("\n------ k = %d ------\n", k);
        printf("beta: %lf\n", beta);
        printf("\n");

        // fork aqui
        for (i=0; i<k; i++) {
            sum += r[k-i-1]*y[i];

            printf("-- i = %d --\n", i);
            printf("sum: %lf\n", sum);
        }
        // join aqui

        alpha = - (r[k] + sum)/beta;
        printf("\nalpha: %lf\n\n", alpha);

        // fork aqui
        for (i=0; i<k; i++) {
            z[i] = y[i] + alpha*y[k-i-1];

            printf("-- i = %d --\n", i);
            printf("z[%d]: %0.2lf\n", i, z[i]);
        }
        printf("\n");

        // barreira aqui

        for (i=0; i<k; i++) {
            y[i] = z[i];
            printf("-- i = %d --\n", i);
            printf("y[%d]: %0.2lf\n", i, y[i]);
        }
        printf("\n");
        // join aqui

        y[k] = alpha;
        printf("y[%d]: %lf\n", k, y[k]);
    }
    #pragma endscop
}


int main(int argc, char** argv) {
    /*Lidar com argumentos passados por CLI*/
    for (int i = 0; i < argc; i++) {
        /* Define tamanho da entrada */
        if (strcmp(argv[i], "-d") == 0) {
            i++;

            if (strcmp(argv[i], "small") == 0) {
                N = 4;  // 250k antes
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

    /* Prevent dead-code elimination. All live-out data must be printed
     *    by the function call in argument. */
    polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

    /* Be clean. */
    POLYBENCH_FREE_ARRAY(r);
    POLYBENCH_FREE_ARRAY(y);

    return 0;
}
