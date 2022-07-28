#include "spmat.h"
#include <stdlib.h> /* malloc, etc */
#include <string.h> /* memcpy */
#include <ctype.h>  /* tolower */
#include <math.h>   /* floor */
#include <limits.h> /* LONG_MAX */
#include <assert.h>

int index_compare(const void *a, const void *b)
{
    index_t ai = *((index_t*)a);
    index_t bi = *((index_t*)b);

    if      (ai < bi) return -1;
    else if (ai > bi) return  1;
    else              return  0;
}

/** spmat_init - Initialize a CSC matrix.
 *
 *  index_t   @m    number of rows
 *  index_t   @n    number of columns
 *  index_t[] @ir   array of row indices (size @jc[@n])
 *  index_t[] @jc   array of column pointers (size @n+1)
 *  num_t[]   @vals array of nonzero values (size @jc[@n]), NULL if pattern matrix
 *
 *  Returns: a newly allocated spmat object.
 *
 *  The returned spmat object takes ownership of all arrays passed to it, so
 *  the user must be careful not to free the arrays he/she allocated for it
 *  except through a call to to spmat_free(). */
spmat *spmat_init(index_t m, index_t n, index_t *ir, index_t *jc, num_t *vals)
{
    spmat *A = malloc(sizeof(spmat));

    A->m    = m;
    A->n    = n;
    A->ir   = ir;
    A->jc   = jc;
    A->vals = vals;

    return A;
}

/** spmat_create - Initialize a CSC matrix using unordered triples.
 *
 *  index_t   @m    number of rows
 *  index_t   @n    number of columns
 *  index_t   @nz   number of nonzeros
 *  index_t[] @ir   array of row pointers (size @nz)
 *  index_t[] @jc   array of column pointers (size @n+1)
 *  num_t[]   @vals array of nonzero values (size @nz), NULL if pattern matrix
 *
 *  Returns: a newly allocated spmat object.
 *
 *  This routine will construct a new CSC matrix using the unordered triples passed
 *  by the user. New arrays will be allocated inside the spmat object, so the
 *  user does not have to worry about freeing the arrays passed to this routine too early. */
spmat *spmat_create(index_t m, index_t n, index_t nz, index_t *ir, index_t *jc, num_t *vals)
{
    index_t k, p, nzcnt, *colcnts;
    spmat *A;

    A = malloc(sizeof(spmat));

    A->m = m;
    A->n = n;
    A->ir = malloc(nz * sizeof(index_t));
    A->jc = calloc(n+1, sizeof(index_t));
    A->vals = vals? malloc(nz * sizeof(num_t)) : NULL;
    colcnts = calloc(n+1, sizeof(index_t));

    for (k = 0; k < nz; ++k)
        ++colcnts[jc[k]];

    nzcnt = 0;

    for (k = 0; k < n; ++k)
    {
        A->jc[k] = nzcnt;
        nzcnt += colcnts[k];
        colcnts[k] = A->jc[k];
    }

    A->jc[n] = nzcnt;

    for (k = 0; k < nz; ++k)
    {
        p = colcnts[jc[k]]++;
        A->ir[p] = ir[k];
        if (vals) A->vals[p] = vals[k];
    }

    free(colcnts);
    return A;
}

/** spmat_copy - Allocates and returns a copy of the passed spmat matrix.
 *
 *  spmat* @A - matrix to copy
 *
 *  Returns: A newly allocated copy of the input matrix.
 */
spmat *spmat_copy(const spmat *A)
{
    index_t m, n, nz, *jc, *ir;
    num_t *vals;

    m    = A->m;
    n    = A->n;
    nz   = spmat_nzs(A);
    jc   = malloc((n+1) * sizeof(index_t));
    ir   = malloc(nz * sizeof(index_t));
    vals = spmat_binary(A)? NULL : malloc(nz * sizeof(num_t));

    memcpy(jc, A->jc, (n+1) * sizeof(index_t));
    memcpy(ir, A->ir, nz * sizeof(index_t));
    if (vals) memcpy(vals, A->vals, nz * sizeof(num_t));

    return spmat_init(m, n, ir, jc, vals);
}

/** spmat_move - Move the contents of one spmat object into another.
 *
 *  spmat* @dest - destination matrix
 *  spmat* @src  - source matrix
 *
 *  This routine assumes that @dest is a spmat object that has been
 *  allocated and owns allocated arrays. */
void spmat_move(spmat *dest, spmat *src)
{
    /* free whatever contents were already in @dest */
    free(dest->ir);
    free(dest->jc);
    free(dest->vals);

    /* copy dimensions and pointers from @src into @dest */
    *dest = *src;

    /* free the @src structure */
    src->m = src->n = 0;
    src->ir = src->jc = NULL, src->vals = NULL;
    free(src);
}

/** spmat_free - Frees memory owned by spmat object and then frees the spmat object itself.
 *
 *  spmat* @A - matrix to free
 */
void spmat_free(spmat *A)
{
    if (!A) return;

    A->m = A->n = 0;

    free(A->jc);
    free(A->ir);
    free(A->vals);
    free(A);
}

spmat* spmat_read(FILE *f)
{
    char line[1025], banner[64], mtx[64], crd[64], data_type[64], storage_scheme[64];

    assert(f && fgets(line, 1025, f));

    assert(sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, storage_scheme) == 5);

    char *p;

    for (p = mtx; *p; *p=tolower(*p), ++p);
    for (p = crd; *p; *p=tolower(*p), ++p);
    for (p = data_type; *p; *p=tolower(*p), ++p);
    for (p = storage_scheme; *p; *p=tolower(*p), ++p);

    assert(!strcmp(banner, "%%MatrixMarket"));
    assert(!strcmp(mtx, "matrix"));
    assert(!strcmp(crd, "coordinate"));

    int binary;

    if (!strcmp(data_type, "real") || !strcmp(data_type, "integer")) binary = 0;
    else if (!strcmp(data_type, "pattern")) binary = 1;
    else assert(0);

    int symmetric;

    if (!strcmp(storage_scheme, "general")) symmetric = 0;
    else if (!strcmp(storage_scheme, "symmetric")) symmetric = 1;
    else assert(0);

    while (fgets(line, 1025, f) && line[0] == '%');

    index_t m, n, nz;

    assert(sscanf(line, "%ld %ld %ld", &m, &n, &nz) == 3);

    if (symmetric) nz *= 2;

    index_t *ir = malloc(nz * sizeof(index_t));
    index_t *jc = malloc(nz * sizeof(index_t));
    num_t *vals = binary? NULL : malloc(nz * sizeof(num_t));

    index_t i, j, k;
    num_t v;

    k = 0;

    while (fgets(line, 1025, f))
    {
        sscanf(line, "%ld %ld %lg", &i, &j, &v);
        ir[k] = i-1;
        jc[k] = j-1;
        if (!binary) vals[k] = v;

        ++k;

        if (symmetric)
        {
            ir[k] = j-1;
            jc[k] = i-1;
            if (!binary) vals[k] = v;
            ++k;
        }
    }

    spmat *A = spmat_create(m, n, nz, ir, jc, vals);

    free(ir);
    free(jc);
    free(vals);

    return A;
}

void spmat_write(const spmat *A, FILE *f)
{
    index_t i, j, p;

    if (spmat_binary(A)) fprintf(f, "%%%%MatrixMarket matrix coordinate pattern general\n");
    else fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");

    fprintf(f, "%ld %ld %ld\n", A->m, A->n, spmat_nzs(A));

    for (j = 0; j < A->n; ++j)
        for (p = A->jc[j]; p < A->jc[j+1]; ++p)
        {
            i = A->ir[p];
            if (!spmat_binary(A)) fprintf(f, "%ld %ld %lg\n", i+1, j+1, A->vals[p]);
            else fprintf(f, "%ld %ld\n", i+1, j+1);
        }
}

spmat *spmat_read_from_file(const char *filename)
{
    FILE *f = fopen(filename, "r");
    spmat *A = spmat_read(f);
    fclose(f);
    return A;
}

void spmat_write_to_file(const spmat *A, const char *filename)
{
    FILE *f = fopen(filename, "w");
    spmat_write(A, f);
    fclose(f);
}

void spmat_pretty(const spmat *A, FILE *f)
{
    spmat *M = spmat_transpose(A);

    char *line = malloc(2*A->n + 1);
    line[2*A->n] = 0;

    index_t i;

    for (index_t j = 0; j < M->n; ++j)
    {
        for (i = 0; i < A->n; ++i)
            line[2*i] = line[2*i+1] = '.';

        for (index_t ip = M->jc[j]; ip < M->jc[j+1]; ++ip)
        {
            i = M->ir[ip];
            line[2*i] = 'X';
        }

        fprintf(f, "%s\n", line);
    }

    free(line);
    spmat_free(M);
}

void spmat_pretty_stdout(const spmat *A)
{
    spmat_pretty(A, stdout);
}

/* don't call this function when (vals != NULL) && (A->vals == NULL) */
void spmat_triples(const spmat *A, index_t **ir, index_t **jc, num_t **vals)
{
    index_t j, p, nz, *_ir, *_jc;
    num_t *_vals;

    nz = spmat_nzs(A);
    _vals = NULL;

    _ir = *ir = malloc(nz * sizeof(index_t));
    _jc = *jc = malloc(nz * sizeof(index_t));
    if (!spmat_binary(A)) _vals = *vals = malloc(nz * sizeof(num_t));

    for (j = 0; j < A->n; ++j)
        for (p = A->jc[j]; p < A->jc[j+1]; ++p)
        {
            *_ir++ = A->ir[p];
            *_jc++ = j;
            if (vals) *_vals++ = A->vals[p];
        }
}

spmat *spmat_transpose(const spmat *A)
{
    spmat *AT;
    index_t *ir, *jc;
    num_t *vals = NULL;

    spmat_triples(A, &jc, &ir, spmat_binary(A)? NULL : &vals);
    AT = spmat_create(A->n, A->m, spmat_nzs(A), ir, jc, vals);

    free(ir);
    free(jc);
    free(vals);

    return AT;
}

void transpose(spmat *A)
{
    spmat *AT = spmat_transpose(A);
    spmat_move(A, AT);
}

typedef enum spa_state_t {SPA_STATE_ALLOWED, SPA_STATE_NOT_ALLOWED, SPA_STATE_SET} spa_state_t;

index_t *spmat_spmv(const spmat *A, index_t *x)
{
    index_t n = A->n;
    index_t m = A->m;

    int binary = spmat_binary(A);

    index_t *y = calloc(m, sizeof(index_t));

    for (index_t j = 0; j < n; ++j)
        if (x[j])
            for (index_t ip = A->jc[j]; ip < A->jc[j+1]; ++ip)
            {
                index_t i = A->ir[ip];
                y[i] = binary? 1 : (y[i] + (x[j] * A->vals[ip]));
            }

    return y;
}

void spmv(const spmat *A, index_t *x)
{
    index_t *y = spmat_spmv(A, x);
    memcpy(x, y, A->m * sizeof(index_t));
    free(y);
}

void bfs_spmv(const spmat *A, index_t *f, index_t *v, index_t *fs)
{
    index_t n = A->n;
    assert(n == A->m);

    #ifdef THREADED
    #pragma omp parallel for
    #endif
    for (index_t i = 0; i < n; ++i)
    {
        fs[i] = 0;

        if (!v[i])
        {
            for (index_t kp = A->jc[i]; kp < A->jc[i+1]; ++kp)
            {
                index_t k = A->ir[kp];

                if (f[k]) fs[i] = 1;
            }
        }
    }
}

spmat *spmat_add(const spmat *A, const spmat *B)
{
    assert(A->m == B->m && A->n == B->n);

    int binaryA = spmat_binary(A);
    int binaryB = spmat_binary(B);
    int binaryC = (binaryA && binaryB);

    index_t m = A->m;
    index_t n = B->n;
    index_t nz = 0;
    index_t nzmax = spmat_nzs(A) + spmat_nzs(B);

    index_t *jc = malloc((n+1) * sizeof(index_t));
    index_t *ir = malloc(nzmax * sizeof(index_t));
    num_t *vals = binaryC? NULL : malloc(nzmax * sizeof(num_t));

    spa_state_t *spa_states = malloc(m * sizeof(spa_state_t));
    num_t *spa_vals = binaryC? NULL : malloc(m * sizeof(num_t));

    for (index_t j = 0; j < n; ++j)
    {
        if (nz + m > nzmax)
        {
            nzmax = nz + 2*m;
            ir = realloc(ir, nzmax * sizeof(index_t));
        }

        jc[j] = nz;

        for (index_t i = 0; i < m; ++i)
            spa_states[i] = SPA_STATE_ALLOWED;

        for (index_t ip = A->jc[j]; ip < A->jc[j+1]; ++ip)
        {
            index_t i = A->ir[ip];
            num_t Aij = binaryA? 1 : A->vals[ip];

            if (!binaryC) spa_vals[i] = Aij;
            spa_states[i] = SPA_STATE_SET;
            ir[nz++] = i;
        }

        for (index_t ip = B->jc[j]; ip < B->jc[j+1]; ++ip)
        {
            index_t i = B->ir[ip];
            num_t Bij = binaryB? 1 : B->vals[ip];

            switch (spa_states[i])
            {
                case SPA_STATE_ALLOWED:
                    if (!binaryC) spa_vals[i] = Bij;
                    spa_states[i] = SPA_STATE_SET;
                    ir[nz++] = i;
                    break;
                case SPA_STATE_SET:
                    if (!binaryC) spa_vals[i] += Bij;
                    break;
                case SPA_STATE_NOT_ALLOWED:
                    assert(0);
            }
        }

        for (index_t p = jc[j]; p < nz; ++p)
        {
            index_t i = ir[p];

            if (!binaryC)
                if (spa_states[i] == SPA_STATE_SET)
                    vals[p] = spa_vals[i];
        }
    }

    jc[n] = nz;

    free(spa_states);
    free(spa_vals);

    ir = realloc(ir, nz * sizeof(index_t));
    if (!binaryC) vals = realloc(vals, nz * sizeof(num_t));

    return spmat_init(m, n, ir, jc, vals);
}

void add(spmat *A, const spmat *B)
{
    spmat *C = spmat_add(A, B);
    spmat_move(A, C);
}

spmat* spmat_spgemm(const spmat *A, const spmat *B, const spmat *M, int notM)
{
    assert(A->n == B->m);

    int binaryA = spmat_binary(A);
    int binaryB = spmat_binary(B);
    int binaryC = (binaryA && binaryB);
    int masked = (M != NULL);

    if (masked) assert(A->m == M->m && B->n == M->n);

    index_t m = A->m;
    index_t n = B->n;
    index_t nz = 0;

    index_t nzmax = masked? (notM? ((M->m * M->n) - spmat_nzs(M)) : spmat_nzs(M)) : spmat_nzs(A) + spmat_nzs(B);

    index_t *jc = malloc((n+1) * sizeof(index_t));
    index_t *ir = malloc(nzmax * sizeof(index_t));
    num_t *vals = binaryC? NULL : malloc(nzmax * sizeof(num_t));

    spa_state_t *spa_states = malloc(m * sizeof(spa_state_t));
    num_t *spa_vals = binaryC? NULL : malloc(m * sizeof(num_t));

    for (index_t j = 0; j < n; ++j)
    {
        jc[j] = nz;

        for (index_t i = 0; i < m; ++i)
            spa_states[i] = SPA_STATE_ALLOWED;

        if (!masked && nz + m > nzmax)
        {
            nzmax = nz + m;
            ir = realloc(ir, nzmax * sizeof(index_t));
            if (!binaryC) vals = realloc(vals, nzmax * sizeof(num_t));
        }

        if (masked)
            for (index_t p = M->jc[j]; p < M->jc[j+1]; ++p)
                spa_states[M->ir[p]] = notM? SPA_STATE_NOT_ALLOWED : SPA_STATE_ALLOWED;

        for (index_t kp = B->jc[j]; kp < B->jc[j+1]; ++kp)
        {
            index_t k = B->ir[kp];
            num_t Bkj = binaryB? 1 : B->vals[kp];

            for (index_t ip = A->jc[k]; ip < A->jc[k+1]; ++ip)
            {
                index_t i = A->ir[ip];
                num_t Aik = binaryA? 1 : A->vals[ip];

                switch (spa_states[i])
                {
                    case SPA_STATE_ALLOWED:
                        if (!binaryC) spa_vals[i] = Aik * Bkj;
                        spa_states[i] = SPA_STATE_SET;
                        ir[nz++] = i;
                        break;
                    case SPA_STATE_SET:
                        if (!binaryC) spa_vals[i] += Aik * Bkj;
                        break;
                    case SPA_STATE_NOT_ALLOWED:
                        break;
                }
            }
        }

        for (index_t p = jc[j]; p < nz; ++p)
        {
            index_t i = ir[p];

            if (!binaryC)
                if (spa_states[i] == SPA_STATE_SET)
                    vals[p] = spa_vals[i];
        }
    }

    free(spa_states);
    free(spa_vals);

    jc[n] = nz;


    ir = realloc(ir, nz * sizeof(index_t));
    if (!binaryC) vals = realloc(vals, nz * sizeof(num_t));

    return spmat_init(m, n, ir, jc, vals);
}

void spgemm(const spmat *A, spmat *B, const spmat *M, int notM)
{
    spmat *C = spmat_spgemm(A, B, M, notM);
    spmat_move(B, C);
}

void ewiseapply(spmat *A, num_t val, const spmat *M)
{
    index_t m = A->m;
    index_t n = A->n;
    index_t nz = spmat_nzs(A);
    int binaryA = spmat_binary(A);

    if (!M)
    {
        if (binaryA) A->vals = malloc(nz * sizeof(num_t));

        for (index_t i = 0; i < nz; ++i)
            A->vals[i] = val;
    }
    else
    {
        assert(A->m == M->m && A->n == M->n);

        spa_state_t *spa_states = malloc(m * sizeof(spa_state_t));

        for (index_t j = 0; j < n; ++j)
        {
            for (index_t p = 0; p < m; ++p)
                spa_states[p] = SPA_STATE_NOT_ALLOWED;

            for (index_t p = M->jc[j]; p < M->jc[j+1]; ++p)
                spa_states[M->ir[p]] = SPA_STATE_ALLOWED;

            for (index_t ip = A->jc[j]; ip < A->jc[j+1]; ++ip)
            {
                index_t i = A->ir[ip];

                if (spa_states[i] == SPA_STATE_ALLOWED)
                    A->vals[ip] = val;
            }
        }

        free(spa_states);
    }
}

index_t dense_vector_nzs(index_t *v, index_t n)
{
    index_t nz = 0;

    #ifdef THREADED
    #pragma omp parallel for reduction(+:nz)
    #endif
    for (index_t i = 0; i < n; ++i)
        if (v[i]) ++nz;

    return nz;
}

void dense_vector_apply(index_t *v, index_t val, const index_t *x, index_t n)
{
    if (!x)
    {
        #ifdef THREADED
        #pragma omp parallel for
        #endif
        for (index_t i = 0; i < n; ++i)
            v[i] = val;
    }
    else
    {
        #ifdef THREADED
        #pragma omp parallel for
        #endif
        for (index_t i = 0; i < n; ++i)
            if (x[i]) v[i] = val;
    }
}

index_t *bfs_fast(const spmat *A, index_t s, index_t *iters)
{
    index_t n = A->m;
    assert(A->m == A->n && s >= 0 && s < n);

    index_t *f  = calloc(n, sizeof(index_t));  // frontier
    index_t *v  = calloc(n, sizeof(index_t));  // visited
    index_t *w  = malloc(n * sizeof(index_t)); // levels
    index_t *fs = malloc(n * sizeof(index_t));  // next frontier

    #ifdef THREADED
    #pragma omp parallel for
    #endif
    for (index_t i = 0; i < n; ++i)
        w[i] = -1;

    index_t l;
    f[s] = 1;
    w[s] = l = 0;

    double t0, t1;
    double t_nzs = 0;
    double t_apply = 0;
    double t_spmv = 0;
    double t_overall = omp_get_wtime();

    for (;;)
    {
        t0 = omp_get_wtime();
        index_t nzs = dense_vector_nzs(f, n);
        t1 = omp_get_wtime();
        t_nzs += (t1-t0);

        if (!nzs) break;

        t0 = omp_get_wtime();
        bfs_spmv(A, f, v, fs); // fs <- ~v && (A^T ||.&& f), v <- v && fs
        t1 = omp_get_wtime();
        t_spmv += (t1-t0);

        index_t *t = f;
        f = fs;
        fs = t;

        t0 = omp_get_wtime();
        dense_vector_apply(v, 1, f, n);
        dense_vector_apply(w, l+1, f, n); // w[f] <- l+1
        t1 = omp_get_wtime();
        t_apply += (t1-t0);

        ++l;
    }

    t_overall = omp_get_wtime() - t_overall;

    int nthreads = 1;

    #ifdef THREADED
    #pragma omp parallel
    {
        if (!omp_get_thread_num()) nthreads = omp_get_num_threads();
    }
    #endif

    tprintf("Performed %ld breadth first search iterations from vertex %ld [%.3f secs] [%d thread(s)]\n", l, s+1, t_overall, nthreads);
    tprintf("Breakdown:\n");
    tprintf("    dense_vector_nzs()   - [%.3f secs]\n", t_nzs);
    tprintf("    dense_vector_apply() - [%.3f secs]\n", t_apply);
    tprintf("    bfs_spmv()           - [%.3f secs]\n\n", t_spmv);

    if (iters) *iters = l;

    free(v);
    free(f);
    free(fs);

    return w;
}

index_t *bfs(const spmat *A, index_t s, index_t *iters)
{
    index_t n = A->m;
    assert(A->m == A->n && s >= 0 && s < n);

    index_t *visited  = calloc(n, sizeof(index_t));
    index_t *frontier = calloc(n, sizeof(index_t));
    index_t *levels   = malloc(n * sizeof(index_t));

    for (index_t i = 0; i < n; ++i)
        levels[i] = -1;

    spmat *AT = spmat_transpose(A);

    index_t level;

    frontier[s] = 1;
    levels[s] = level = 0;

    double t0, t1;
    double t_nzs = 0;
    double t_apply = 0;
    double t_spmv = 0;
    double t_overall = omp_get_wtime();

    for (;;)
    {
        t0 = omp_get_wtime();
        index_t nzs = dense_vector_nzs(frontier, n);
        t1 = omp_get_wtime();
        t_nzs += (t1-t0);

        if (!nzs) break;

        t0 = omp_get_wtime();
        dense_vector_apply(visited, 1, frontier, n);
        t1 = omp_get_wtime();
        t_apply += (t1-t0);

        t0 = omp_get_wtime();
        spmv(AT, frontier);
        t1 = omp_get_wtime();
        t_spmv += (t1-t0);

        t0 = omp_get_wtime();
        dense_vector_apply(frontier, 0, visited, n);
        dense_vector_apply(levels, level+1, frontier, n);
        t1 = omp_get_wtime();
        t_apply += (t1-t0);

        printf("%ld\n", level);

        ++level;
    }

    t_overall = omp_get_wtime() - t_overall;

    int nthreads = 1;

    #ifdef THREADED
    #pragma omp parallel
    {
        if (!omp_get_thread_num()) nthreads = omp_get_num_threads();
    }
    #endif

    tprintf("Performed %ld breadth first search iterations from vertex %ld [%.3f secs] [%d thread(s)]\n", level, s+1, t_overall, nthreads);
    tprintf("Breakdown:\n");
    tprintf("    dense_vector_nzs()   - [%.3f secs]\n", t_nzs);
    tprintf("    dense_vector_apply() - [%.3f secs]\n", t_apply);
    tprintf("    spmv()               - [%.3f secs]\n\n", t_spmv);

    if (iters) *iters = level;

    free(visited);
    free(frontier);
    spmat_free(AT);

    return levels;
}

index_t *ullman_yannakakis(const spmat *A, index_t s, index_t *iters)
{
    index_t n = A->m;
    assert(A->m == A->n && s >= 0 && s < n);

    double t0, t1;

    t0 = omp_get_wtime();
    index_t k = (index_t)floor(2*sqrt(n)); /* k = O(n^(1/2)log(n)) is the number of distinguished vertices */
    index_t *d = malloc((k+1)*sizeof(index_t));

    /* randomly sample k vertices from [0..n-1] without replacement, using Fisher-Yates method */
    index_t *x = malloc(n * (sizeof(index_t)));

    index_t l;

    for (l = 0; l < n; ++l)
        x[l] = l;

    for (index_t i = 0; i < k; ++i)
    {
        index_t r = rand() % (n-i);
        index_t t = x[r];
        x[r] = x[n-i-1];
        x[n-i-1] = t;
    }

    int contains_source = 0;

    for (index_t i = 0; i < k; ++i)
        if ((d[i] = x[n-k+i]) == s)
            contains_source = 1; /* distinguished set already includes source */

    free(x);

    if (!contains_source)
        d[k++] = s; /* add source to distinguished set */

    t1 = omp_get_wtime();
    vtprintf("\nRandomly sampled %ld distinguished vertices [%.3f secs]\n", k, t1-t0);

    t0 = omp_get_wtime();
    qsort(d, k, sizeof(index_t), index_compare);
    t1 = omp_get_wtime();
    vtprintf("\nSorted %ld integers                         [%.3f secs]\n", k, t1-t0);

    t0 = omp_get_wtime();
    index_t *ir = malloc(k * sizeof(index_t));
    index_t *jc = malloc((k+1) * sizeof(index_t));

    for (index_t i = 0; i < k; ++i)
    {
        jc[i] = i;
        ir[i] = d[i];
    }

    jc[k] = k;

    spmat *F = spmat_init(n, k, ir, jc, NULL); /* n by sqrt(n)log(n) frontier matrix */
    spmat *P = spmat_copy(F); /* bfs levels matrix */

    ewiseapply(P, 0, NULL);

    spmat *AT = spmat_transpose(A);

    spmat *N;

    index_t bfsiters = (index_t)floor(sqrt(n)*log(n)); /* limited search count */

    t1 = omp_get_wtime();
    tprintf("\nInitialized frontier, visited, and levels matrices [%.3f secs]\n\n", t1-t0);

    t0 = omp_get_wtime();
    for (index_t i = 0; i < bfsiters; ++i)
    {
        t1 = omp_get_wtime();
        N = spmat_spgemm(AT, F, P, 1);

        vtprintf(" 1/4. N <- V .* (A^T * F), %ld/%ld iters performed [%.3f secs]\n", i, bfsiters, omp_get_wtime()-t1);

        add(P, N);
        vtprintf(" 2/4. P <- (P + N),        %ld/%ld iters performed [%.3f secs]\n", i, bfsiters, omp_get_wtime()-t1);

        ewiseapply(P, (num_t)i+1, N);
        vtprintf(" 3/4. P<N> <- i+1,         %ld/%ld iters performed [%.3f secs]\n", i, bfsiters, omp_get_wtime()-t1);

        spmat_move(F, N);
        double t2 = omp_get_wtime();

        vtprintf(" 4/4. F <- N,              %ld/%ld iters performed [%.3f secs]\n\n", i, bfsiters, t2-t1);
        tprintf("%ld/%ld iters performed [%.3f secs]\n\n", i+1, bfsiters, t2-t0);
    }
    t1 = omp_get_wtime();

    tprintf("Performed %ld breadth first search iterations from each distinguished vertex [%f secs]\n", bfsiters, t1-t0);

    index_t *mapback = calloc(n, sizeof(index_t)); /* 0 means not in distinguished set, nonzero i means i-1 is [0..k-1] index */

    for (index_t i = 0; i < k; ++i)
        mapback[d[i]] = i+1;

    index_t *D = malloc((k*k) * sizeof(num_t)); /* all pairs shortest paths distance matrix for distinguished vertices */

    for (index_t i = 0; i < k*k; ++i)
        D[i] = -1;

    for (index_t j = 0; j < k; ++j)
    {
        for (index_t ip = P->jc[j]; ip < P->jc[j+1]; ++ip)
        {
            index_t i = P->ir[ip];
            index_t distance = (index_t)P->vals[ip];

            if (mapback[i])
                D[(mapback[i]-1)*k + j] = distance;
        }
    }

    /* All-pairs shortest path (Floyd-Warshall) */
    for (index_t p = 0; p < k; ++p)
        for (index_t i = 0; i < k; ++i)
            for (index_t j = 0; j < k; ++j)
                if (i != j)
                {
                    index_t Dij = D[i*k + j];
                    index_t Dip = D[i*k + p];
                    index_t Dpj = D[p*k + j];

                    if (Dip != -1 && Dpj != -1)
                        D[i*k + j] = MIN(Dij, Dip + Dpj);
                }

    index_t *levels = malloc(n * sizeof(index_t));

    for (index_t i = 0; i < n; ++i)
        levels[i] = -1;

    /* computing estimated shortest path from s to v (s = source) */

    for (index_t v = 0; v < n; ++v)
    {
        index_t stv = LONG_MAX;

        /* From Klein paper:
         *
         * "Take minimum of all distinguished vertices x, of the
         *  auxiliary-graph shortest path length from s to x plus
         *  the minimum path-size from x to v found in x's limited
         *  search." */
        for (index_t x = 0; x < k; ++x)
        {
            /* get auxiliary-graph shortest path length from s to x */
            index_t stx = D[(mapback[s]-1)*k + x];

            index_t xtv = -1;
            /* get minimum path-size from x to v found in x's limited
             * search */
            for (index_t ip = P->jc[x]; ip < P->jc[x+1]; ++ip)
            {
                index_t i = P->ir[ip];
                index_t distance = (index_t)P->vals[ip];

                if (i == v)
                {
                    xtv = distance;
                    break;
                }
            }

            if (stx != -1 && xtv != -1)
                stv = MIN(stv, stx + xtv);
        }

        levels[v] = (stv < LONG_MAX)? stv : -1;
    }

    *iters = bfsiters;

    free(D);
    free(mapback);
    free(d);

    return levels;

}
