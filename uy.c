#include "spmat.h"
#include <stdio.h>
#include <stdlib.h> /* strtol */
#include <time.h>
#include <unistd.h>

double t0, t1;

int main(int argc, char *argv[])
{
    srand(time(0)*getpid());

    if (argc < 4)
    {
        fprintf(stderr, "usage: %s <mtx> <paths> [sources...]\n", *argv);
        return -1;
    }

    const char *mtx = argv[1];
    const char *paths = argv[2];

    tprintf("Started reading '%s'\n", mtx);
    t0 = omp_get_wtime();
    spmat *A = spmat_read_from_file(mtx);
    t1 = omp_get_wtime();

    tprintf("Finished reading '%s' (%ld rows, %ld columns, %ld nonzeros) [%.3f secs]\n\n", mtx, A->m, A->n, spmat_nzs(A), t1-t0);

    FILE *f = fopen(paths, "w");

    for (int argidx = 3; argidx < argc; ++argidx)
    {
        index_t iters;
        index_t s = strtol(argv[argidx], NULL, 10) - 1;

        tprintf("Starting ullman-yannakakis from vertex %ld\n", s+1);
        index_t *levels = ullman_yannakakis(A, s, &iters);

        #ifdef LOGGER
        index_t explored = 0;
        #pragma omp parallel for reduction(+:explored)
        for (index_t i = 0; i < A->n; ++i)
            if (levels[i] != -1) ++explored;
        fprintf(stderr, "(at least) %ld/%ld vertices reachable from vertex %ld\n\n", explored, A->n, s+1);
        #endif

        fprintf(f, "** ullman-yannakakis levels from vertex %ld **\n", s+1);
        for (index_t i = 0; i < A->n-1; ++i)
            fprintf(f, "%ld ", levels[i]);
        fprintf(f, "%ld\n\n", levels[A->n-1]);

        free(levels);
    }

    fclose(f);
    spmat_free(A);
    return 0;
}
