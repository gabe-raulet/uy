#ifndef ULLMAN_YANNAKAKIS_H_
#define ULLMAN_YANNAKAKIS_H_

#include <stdio.h>

#define MAX(a, b) (((a) < (b))? (b) : (a))
#define MIN(a, b) (((a) < (b))? (a) : (b))

#ifdef LOGGER
#include <omp.h>
#define timer_start() do { t0 = omp_get_wtime(); } while (0)
#define timer_stop() do { t1 = omp_get_wtime(); } while (0)
#define timer_elapsed() (t1-t0+0.0)
#define tprintf(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)
#else
#define timer_start()
#define timer_stop()
#define timer_elapsed()
#define tprintf(fmt, ...)
#endif

#ifdef THREADED
#include <omp.h>
#endif

typedef signed long index_t;
typedef double num_t;

typedef struct spmat
{
    index_t m;   /* number of rows */
    index_t n;   /* number of columns */
    index_t *ir; /* row indices */
    index_t *jc; /* column pointers */
    num_t *vals; /* nonzeros (NULL if pattern matrix) */
} spmat;

#define spmat_nzs(A) ((A)->jc[(A)->n])
#define spmat_binary(A) ((A)->vals==NULL)

/* workhorse routines */
spmat*   spmat_init        (index_t m, index_t n, index_t *ir, index_t *jc, num_t *vals);
spmat*   spmat_create      (index_t m, index_t n, index_t nz, index_t *ir, index_t *jc, num_t *vals);
spmat*   spmat_copy        (const spmat *A);
void     spmat_move        (spmat *dest, spmat *src);
void     spmat_free        (spmat *A);
spmat*   spmat_read        (FILE *f);
void     spmat_write       (const spmat *A, FILE *f);
void     spmat_triples     (const spmat *A, index_t **ir, index_t **jc, num_t **vals);
spmat*   spmat_transpose   (const spmat *A);
index_t* spmat_spmv        (const spmat *A, index_t *x);
spmat*   spmat_add         (const spmat *A, const spmat *B);
spmat*   spmat_spgemm      (const spmat *A, const spmat *B, const spmat *M, int notM);
index_t* bfs               (const spmat *A, index_t s, index_t *iters);
index_t* ullman_yannakakis (const spmat *A, index_t s, index_t *iters);

/* input/output subroutines */
spmat* spmat_read_from_file (const char *filename);
void   spmat_write_to_file  (const spmat *A, const char *filename);
void   spmat_pretty         (const spmat *A, FILE *f);
void   spmat_pretty_stdout  (const spmat *A);

/* algebraic subroutines */
void transpose  (spmat *A);                                           /* A <- A^T */
void spmv       (const spmat *A, index_t *x);                         /* x <- A*x */
void add        (spmat *A, const spmat *B);                           /* A <- A+B */
void spgemm     (const spmat *A, spmat *B, const spmat *M, int notM); /* B <- M .* (A*B) */
void ewiseapply (spmat *A, num_t val, const spmat *M);

/* for sorting index_t arrays */
int index_compare(const void *a, const void *b);

#endif
