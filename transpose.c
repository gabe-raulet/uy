#include "uy.h"
#include <unistd.h> /* getopt */

const char *usage =
    "transpose -- Read a sparse matrix and output its transpose\n"
    "Flags:\n"
    " -i -- input matrix filename  (stdin)\n"
    " -o -- output matrix filename (stdout)\n"
    " -h -- help\n";

double t0, t1;

int main(int argc, char *argv[])
{
    char *ifname = NULL;
    char *ofname = NULL;

    extern char *optarg;
    const char *optstring = "hi:o:";
    int c;

    while ((c = getopt(argc, argv, optstring)) != -1)
    {
        switch (c)
        {
            case 'h':
                fprintf(stderr, "%s", usage);
                return -1;
            case 'i': ifname = optarg; break;
            case 'o': ofname = optarg; break;
        }
    }

    FILE *fin = stdin;
    FILE *fout = stdout;

    int fi = 0;
    int fo = 0;

    if (ifname) fin = fopen(ifname, "r"), fi = 1;
    if (ofname) fout = fopen(ofname, "w"), fo = 1;

    tprintf("Started reading '%s'\n", fi? ifname : "stdin");
    timer_start();
    spmat *A = spmat_read(fin);
    timer_stop();
    tprintf("Finished reading '%s' (%ld rows, %ld columns, %ld nonzeros) [%f secs]\n", fi? ifname : "stdin", A->m, A->n, spmat_nzs(A), timer_elapsed());

    tprintf("Computing transpose\n");
    timer_start();
    transpose(A);
    timer_stop();
    tprintf("Computed transpose [%f secs]\n", timer_elapsed());

    spmat_write(A, fout);
    spmat_free(A);

    if (fi) fclose(fin);
    if (fo) fclose(fout);

    return 0;
}
