#include "spmat.h"
#include <stdio.h>
#include <stdlib.h> /* strtol */

#ifdef LOGGER
double t0, t1;
#endif

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        fprintf(stderr, "usage: %s <mtx> <paths> [sources...]\n", *argv);
        return -1;
    }

    const char *mtx = argv[1];
    const char *paths = argv[2];

    tprintf("Started reading '%s'\n", mtx);
    timer_start();
    spmat *A = spmat_read_from_file(mtx);
    timer_stop();
    tprintf("Finished reading '%s' (%ld rows, %ld columns, %ld nonzeros) [%f secs]\n\n", mtx, A->m, A->n, spmat_nzs(A), timer_elapsed());

    FILE *f = fopen(paths, "w");

    for (int argidx = 3; argidx < argc; ++argidx)
    {
        index_t iters;
        index_t s = strtol(argv[argidx], NULL, 10) - 1;

        tprintf("Starting bfs from vertex %ld\n", s+1);
        timer_start();
        index_t *levels = bfs(A, s, &iters);
        timer_stop();
        tprintf("Performed %ld bfs iterations [%f secs]\n", iters, timer_elapsed());

        #ifdef LOGGER
        index_t explored = 0;
        #pragma omp parallel for reduction(+:explored)
        for (index_t i = 0; i < A->n; ++i)
            if (levels[i] != -1) ++explored;
        fprintf(stderr, "%ld/%ld vertices reachable from vertex %ld\n\n", explored, A->n, s+1);
        #endif

        fprintf(f, "** bfs levels from vertex %ld (eccentricity = %ld) **\n", s+1, iters);
        for (index_t i = 0; i < A->n-1; ++i)
            fprintf(f, "%ld ", levels[i]);
        fprintf(f, "%ld\n\n", levels[A->n-1]);

        free(levels);
    }

    fclose(f);
    spmat_free(A);
    return 0;
}
