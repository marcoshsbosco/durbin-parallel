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
#include <pthread.h>

/* Include polybench common header. */
#include "polybench.h"

/* Include benchmark-specific header. */
#include "durbin.h"


/* Variáveis que tínhamos no arquivo header antes de passar o tamanho do dataset como argumento de CLI */
int N;
int _PB_N;

pthread_barrier_t barrier;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
DATA_TYPE somas_parciais = 0;

struct ArgsFluxo {  // encapsula argumentos para cada fluxo
    int fluxos;
    int fluxoAtual;

    int n;
    DATA_TYPE *r, *y;
};

/* Mensagem de ajuda do programa (--help) */
void help_message() {
    printf("--- Instruções de uso ---\n");
    printf("./durbin.out [OPÇÕES]\n");
    printf("\n[OPÇÕES]\n");
    printf("-d: Tamanho do dataset. Pode ser um dentre {small, medium, large}.\n");
    printf("-t: Quantidade de threads/fluxos.\n");
}

/* Array initialization. */
static void init_array(int n, DATA_TYPE *r) {
    int i, j;

    for (i = 0; i < n; i++) {
        r[i] = (n+1-i);
    }
}

/* DCE code. Must scan the entire live-out data.
 *Can be used also to check the correctness of the output. */
static void print_array(int n, DATA_TYPE *y) {
    int i;

    for (i = 0; i < n; i++) {
        printf(DATA_PRINTF_MODIFIER, y[i]);
    }

    printf("\n");
}


/* Main computational kernel. The whole function will be timed,
 *including the call and return. */
static void *kernel_durbin(void *arg) {
    struct ArgsFluxo args = *(struct ArgsFluxo *)arg;
    DATA_TYPE z[N];
    DATA_TYPE alpha;
    DATA_TYPE beta;
    DATA_TYPE sum;
    int i, k;

    #pragma scop
    args.y[0] = -args.r[0];
    beta = SCALAR_VAL(1.0);
    alpha = -args.r[0];

    for (k = 1; k < args.n; k++) {
        beta = (1-alpha*alpha)*beta;
        sum = SCALAR_VAL(0.0);
        somas_parciais = 0;

        pthread_barrier_wait(&barrier);

        for (i = args.fluxoAtual * k / args.fluxos; i < (k + k * args.fluxoAtual) / args.fluxos; i++) {
            sum += args.r[k-i-1]*args.y[i];
        }

        pthread_mutex_lock(&mutex);
        somas_parciais += sum;
        pthread_mutex_unlock(&mutex);

        pthread_barrier_wait(&barrier);

        sum = somas_parciais;

        alpha = - (args.r[k] + sum)/beta;

        for (i = args.fluxoAtual * k / args.fluxos; i < (k + k * args.fluxoAtual) / args.fluxos; i++) {
            z[i] = args.y[i] + alpha*args.y[k-i-1];
        }

        pthread_barrier_wait(&barrier);

        for (i = args.fluxoAtual * k / args.fluxos; i < (k + k * args.fluxoAtual) / args.fluxos; i++) {
            args.y[i] = z[i];
        }

        args.y[k] = alpha;
    }
    #pragma endscop
}


int main(int argc, char** argv) {
    pthread_t *idsFluxo;
    struct ArgsFluxo *args;
    int fluxos;

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
        } else if (strcmp(argv[i], "-t") == 0) {  // Quantidade de fluxos
            i++;

            fluxos = strtol(argv[i], NULL, 10);
        }
    }

    pthread_barrier_init(&barrier, NULL, fluxos);

    // alocação de memória
    idsFluxo = calloc(fluxos, sizeof(pthread_t));
    args = calloc(fluxos, sizeof(struct ArgsFluxo));

    /* Retrieve problem size. */
    int n = N;
    _PB_N = POLYBENCH_LOOP_BOUND(N,n);

    /* Variable declaration/allocation. */
    DATA_TYPE r[n], y[n];

    /* Initialize array(s). */
    init_array (n, r);

    /* Start timer. */
    polybench_start_instruments;

    // criação de fluxos e de seus atributos (ID e número)
    for (int i = 0; i < fluxos; i++) {
        args[i].fluxoAtual = i;
        args[i].fluxos = fluxos;

        args[i].n = n;

        args[i].r = r;
        args[i].y = y;

        pthread_create(&idsFluxo[i], NULL, &kernel_durbin, &args[i]);
    }

    // join
    for (int i = 0; i < fluxos; i++) {
        pthread_join(idsFluxo[i], NULL);
    }

    pthread_barrier_destroy(&barrier);

    /* Stop and print timer. */
    polybench_stop_instruments;
    polybench_print_instruments;

    return 0;
}