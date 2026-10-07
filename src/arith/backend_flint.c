/**
 * \file   backend_flint.c
 * \brief  FLINT fq_default backend for the pluggable arithmetic interface.
 *
 * Every ring element is an element of R_q = (Z/qZ)[x]/(x^d + 1) held as an
 * fq_default_mat entry; vectors are n x 1 matrices and plain matrices are
 * r x c matrices over the same field.  Ring construction from (q,d) alone:
 * the optional `spec` is ignored.
 */

#include "arith.h"

#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>

#include <stdlib.h>

struct arith_ring
{
  fq_default_ctx_struct *fq;
  fmpz_mod_ctx_struct *mod;
  unsigned long long q;
  unsigned int d;
};

struct arith_poly
{
  struct arith_ring *r;
  fq_default_mat_struct *m; /* 1 x 1 */
};

struct arith_vec
{
  struct arith_ring *r;
  fq_default_mat_struct *m; /* n x 1 */
};

struct arith_mat
{
  struct arith_ring *r;
  fq_default_mat_struct *m; /* nrows x ncols */
};

static inline fq_default_ctx_struct *
ctx_of (const struct arith_ring *r)
{
  return r->fq;
}

static inline fmpz_mod_ctx_struct *
mod_of (const struct arith_ring *r)
{
  return r->mod;
}

static fq_default_mat_struct *
new_mat (struct arith_ring *r, slong rows, slong cols)
{
  fq_default_mat_struct *m = flint_malloc (sizeof (fq_default_mat_struct));
  fq_default_mat_init (m, rows, cols, ctx_of (r));
  return m;
}

static void
free_mat (struct arith_ring *r, fq_default_mat_struct *m)
{
  fq_default_mat_clear (m, ctx_of (r));
  flint_free (m);
}

static void
mat_set_coeff (struct arith_ring *r, fq_default_mat_struct *m, slong row,
               slong col, unsigned int k, long long c)
{
  fq_default_ctx_struct *fq = ctx_of (r);
  fmpz_mod_ctx_struct *mod = mod_of (r);

  fmpz_mod_poly_t p;
  fmpz_mod_poly_init (p, mod);

  fq_default_t e;
  fq_default_init2 (e, fq);

  fq_default_mat_entry (e, m, row, col, fq);
  fq_default_get_fmpz_mod_poly (p, e, fq);
  fmpz_mod_poly_set_coeff_si (p, k, c, mod);
  fq_default_set_fmpz_mod_poly (e, p, fq);
  fq_default_mat_entry_set (m, row, col, e, fq);

  fq_default_clear (e, fq);
  fmpz_mod_poly_clear (p, mod);
}

static long long
mat_get_coeff (struct arith_ring *r, const fq_default_mat_struct *m, slong row,
               slong col, unsigned int k)
{
  fq_default_ctx_struct *fq = ctx_of (r);
  fmpz_mod_ctx_struct *mod = mod_of (r);

  long long c = 0;

  fmpz_mod_poly_t p;
  fmpz_mod_poly_init (p, mod);

  fq_default_t e;
  fq_default_init2 (e, fq);

  fq_default_mat_entry (e, m, row, col, fq);
  fq_default_get_fmpz_mod_poly (p, e, fq);

  fmpz_t v;
  fmpz_init (v);
  fmpz_mod_poly_get_coeff_fmpz (v, p, k, mod);
  c = (long long)fmpz_get_si (v);
  fmpz_clear (v);

  fq_default_clear (e, fq);
  fmpz_mod_poly_clear (p, mod);

  return c;
}

/* ------------------------------------------------------------------ */
/* ring                                                                */
/* ------------------------------------------------------------------ */

static arith_ring_ptr
fl_ring_new (unsigned long long q, unsigned int d, const void *spec)
{
  (void)spec;

  fmpz_t qf;
  fmpz_init_set_ui (qf, q);

  fmpz_mod_ctx_struct *mod = flint_malloc (sizeof (fmpz_mod_ctx_struct));
  fmpz_mod_ctx_init (mod, qf);

  fmpz_mod_poly_t xdplus1;
  fmpz_mod_poly_init (xdplus1, mod);
  fmpz_mod_poly_set_coeff_ui (xdplus1, d, 1, mod);
  fmpz_mod_poly_set_coeff_ui (xdplus1, 0, 1, mod);

  fq_default_ctx_struct *fq = flint_malloc (sizeof (fq_default_ctx_struct));
  fq_default_ctx_init_modulus (fq, xdplus1,
                               mod, "x");

  fmpz_mod_poly_clear (xdplus1, mod);
  fmpz_clear (qf);

  arith_ring_ptr r = (arith_ring_ptr)malloc (sizeof (struct arith_ring));
  r->fq = fq;
  r->mod = mod;
  r->q = q;
  r->d = d;
  return r;
}

static void
fl_ring_free (arith_ring_ptr r)
{
  fq_default_ctx_clear (r->fq);
  fmpz_mod_ctx_clear (r->mod);
  flint_free (r->fq);
  flint_free (r->mod);
  free (r);
}

static unsigned long long
fl_ring_q (const arith_ring_ptr r)
{
  return r->q;
}

static unsigned int
fl_ring_d (const arith_ring_ptr r)
{
  return r->d;
}

/* ------------------------------------------------------------------ */
/* poly                                                                */
/* ------------------------------------------------------------------ */

static arith_poly_ptr
fl_poly_alloc (const arith_ring_ptr r)
{
  arith_poly_ptr p = (arith_poly_ptr)malloc (sizeof (struct arith_poly));
  p->r = r;
  p->m = new_mat (r, 1, 1);
  return p;
}

static void
fl_poly_free (arith_poly_ptr p)
{
  free_mat (p->r, p->m);
  free (p);
}

static void
fl_poly_set_coeff (arith_poly_ptr p, unsigned int i, long long c)
{
  mat_set_coeff (p->r, p->m, 0, 0, i, c);
}

static long long
fl_poly_get_coeff (const arith_poly_ptr p, unsigned int i)
{
  return mat_get_coeff (p->r, p->m, 0, 0, i);
}

/* ------------------------------------------------------------------ */
/* vec                                                                */
/* ------------------------------------------------------------------ */

static arith_vec_ptr
fl_vec_alloc (const arith_ring_ptr r, unsigned int n)
{
  arith_vec_ptr v = (arith_vec_ptr)malloc (sizeof (struct arith_vec));
  v->r = r;
  v->m = new_mat (r, n, 1);
  return v;
}

static void
fl_vec_free (arith_vec_ptr v)
{
  free_mat (v->r, v->m);
  free (v);
}

static void
fl_vec_set (arith_vec_ptr dst, const arith_vec_ptr src)
{
  fq_default_mat_set (dst->m, src->m,
                      ctx_of (dst->r));
}

static void
fl_vec_set_coeff (arith_vec_ptr v, unsigned int elem, unsigned int k,
                  long long c)
{
  mat_set_coeff (v->r, v->m, elem, 0, k, c);
}

static long long
fl_vec_get_coeff (const arith_vec_ptr v, unsigned int elem, unsigned int k)
{
  return mat_get_coeff (v->r, v->m, elem, 0, k);
}

/* ------------------------------------------------------------------ */
/* mat                                                                */
/* ------------------------------------------------------------------ */

static arith_mat_ptr
fl_mat_alloc (const arith_ring_ptr r, unsigned int nrows, unsigned int ncols)
{
  arith_mat_ptr A = (arith_mat_ptr)malloc (sizeof (struct arith_mat));
  A->r = r;
  A->m = new_mat (r, nrows, ncols);
  return A;
}

static void
fl_mat_free (arith_mat_ptr A)
{
  free_mat (A->r, A->m);
  free (A);
}

static void
fl_mat_set_coeff (arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k, long long c)
{
  mat_set_coeff (A->r, A->m, i, j, k, c);
}

static long long
fl_mat_get_coeff (const arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k)
{
  return mat_get_coeff (A->r, A->m, i, j, k);
}

/* ------------------------------------------------------------------ */
/* algebra                                                            */
/* ------------------------------------------------------------------ */

static void
fl_vec_addmul (const arith_vec_ptr w, const arith_mat_ptr A,
               const arith_vec_ptr v)
{
  fq_default_ctx_struct *fq = ctx_of (w->r);

  slong nrows = fq_default_mat_nrows (w->m, fq);
  fq_default_mat_struct *tmp = flint_malloc (sizeof (fq_default_mat_struct));
  fq_default_mat_init (tmp, nrows, 1, fq);

  fq_default_mat_mul (tmp, A->m,
                      v->m, fq);
  fq_default_mat_add (w->m, w->m,
                      tmp, fq);

  fq_default_mat_clear (tmp, fq);
  flint_free (tmp);
}

static void
fl_vec_acc_polyscalar (const arith_vec_ptr w, const arith_poly_ptr s,
                       const arith_vec_ptr v)
{
  fq_default_ctx_struct *fq = ctx_of (w->r);

  fq_default_t e;
  fq_default_init2 (e, fq);
  fq_default_mat_entry (e, s->m, 0, 0, fq);

  fq_default_mat_struct *tmp = flint_malloc (sizeof (fq_default_mat_struct));
  fq_default_mat_init (tmp,
                       fq_default_mat_nrows (v->m, fq),
                       fq_default_mat_ncols (v->m, fq), fq);

  fq_default_mat_scalar_mul (tmp, v->m, e,
                             fq);
  fq_default_mat_add (w->m, w->m,
                      tmp, fq);

  fq_default_mat_clear (tmp, fq);
  flint_free (tmp);
  fq_default_clear (e, fq);
}

static void
fl_vec_acc (const arith_vec_ptr w, const arith_vec_ptr a, int neg)
{
  fq_default_ctx_struct *fq = ctx_of (w->r);

  if (neg)
    fq_default_mat_sub (w->m, w->m,
                        a->m, fq);
  else
    fq_default_mat_add (w->m, w->m,
                        a->m, fq);
}

/* ------------------------------------------------------------------ */

static const arith_ops_s ops = {
    "flint",

    fl_ring_new,
    fl_ring_free,
    fl_ring_q,
    fl_ring_d,

    fl_poly_alloc,
    fl_poly_free,
    fl_poly_set_coeff,
    fl_poly_get_coeff,

    fl_vec_alloc,
    fl_vec_free,
    fl_vec_set,
    fl_vec_set_coeff,
    fl_vec_get_coeff,

    fl_mat_alloc,
    fl_mat_free,
    fl_mat_set_coeff,
    fl_mat_get_coeff,

    fl_vec_addmul,
    fl_vec_acc_polyscalar,
    fl_vec_acc,
};

arith_ops_ptr
arith_flint_ops (void)
{
  return &ops;
}