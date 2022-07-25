#include "spmat.h"
#include <stdio.h>
#include <stdlib.h> /* strtol */

double t0, t1;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        fprintf(stderr,"usage: %s <mtxA> <mtxB> <mtxC = mtxA * mtxB>\n", *argv);
        return -1;
    }

    const char *mtxA = argv[1];
    const char *mtxB = argv[2];
    const char *mtxC = argv[3];

    tprintf("Started reading '%s'\n", mtxA);
    t0 = omp_get_wtime();
    spmat *A = spmat_read_from_file(mtxA);
    t1 = omp_get_wtime();
    tprintf("Finished reading '%s' (%ld rows, %ld columns, %ld nonzeros) [%.3f secs]\n\n", mtxA, A->m, A->n, spmat_nzs(A), t1-t0);

    tprintf("Started reading '%s'\n", mtxA);
    t0 = omp_get_wtime();
    spmat *B = spmat_read_from_file(mtxA);
    t1 = omp_get_wtime();
    tprintf("Finished reading '%s' (%ld rows, %ld columns, %ld nonzeros) [%.3f secs]\n\n", mtxB, B->m, B->n, spmat_nzs(B), t1-t0);

    spmat *C;
    int threads[] = {1, 2, 3, 4, 6, 12};

    for (int i = 0; i < 6; ++i)
    {
        omp_set_num_threads(threads[i]);
        t0 = omp_get_wtime();
        C = spmat_spgemm(A, B, NULL, 0);
        t1 = omp_get_wtime();
        tprintf("C = A * B (%ld rows, %ld columns, %ld nonzeros) [%.3f secs] [%d threads]\n", C->m, C->n, spmat_nzs(C), t1-t0, threads[i]);

        if (i == 6) spmat_write_to_file(C, mtxC);

        spmat_free(C);
    }

    return 0;
}
