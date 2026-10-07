/**
 * \file   backend_lazer.c
 * \brief  Native lazer backend for the pluggable arithmetic interface.
 *
 * \note   The lazer ring object is code-generated/static and can never be
 *         built from (q,d) at runtime.  The optional `spec` handed to
 *         ring_new must therefore point at a compiled-in polyring instance,
 *         otherwise ring_new fails.  All other construction (vec/mat/poly)
 *         is derived from that ring.
 */

#include "arith.h"
#include "../../lazer.h"

#include <stdlib.h>

struct arith_ring
{
  polyring_srcptr ring;
  unsigned long long q;
  unsigned int d;
};

struct arith_poly
{
  arith_ring_ptr ring;
  poly_t p;
};

struct arith_vec
{
  arith_ring_ptr ring;
  polyvec_t v;
};

struct arith_mat
{
  arith_ring_ptr ring;
  polymat_t A;
};

/* ------------------------------------------------------------------ */
/* ring                                                                */
/* ------------------------------------------------------------------ */

static arith_ring_ptr
lz_ring_new (unsigned long long q, unsigned int d, const void *spec)
{
  if (spec == NULL)
    return NULL;

  polyring_srcptr ring = (polyring_srcptr)spec;

  if ((unsigned long long)int_get_i64 (polyring_get_mod (ring)) != q
      || polyring_get_deg (ring) != d)
    return NULL;

  arith_ring_ptr r = (arith_ring_ptr)malloc (sizeof (struct arith_ring));
  r->ring = ring;
  r->q = q;
  r->d = d;
  return r;
}

static void
lz_ring_free (arith_ring_ptr r)
{
  free (r);
}

static unsigned long long
lz_ring_q (const arith_ring_ptr r)
{
  return r->q;
}

static unsigned int
lz_ring_d (const arith_ring_ptr r)
{
  return r->d;
}

/* ------------------------------------------------------------------ */
/* poly                                                               */
/* ------------------------------------------------------------------ */

static arith_poly_ptr
lz_poly_alloc (const arith_ring_ptr r)
{
  arith_poly_ptr p = (arith_poly_ptr)malloc (sizeof (struct arith_poly));
  p->ring = r;
  poly_alloc (p->p, r->ring);
  poly_set_zero (p->p);
  return p;
}

static void
lz_poly_free (arith_poly_ptr p)
{
  poly_free (p->p);
  free (p);
}

static void
lz_poly_set_coeff (arith_poly_ptr p, unsigned int i, long long c)
{
  int_ptr coeff = poly_get_coeff (p->p, i);
  int_set_i64 (coeff, c);
}

static long long
lz_poly_get_coeff (const arith_poly_ptr p, unsigned int i)
{
  long long c = int_get_i64 (poly_get_coeff (p->p, i));
  long long q = (long long)p->ring->q;
  c %= q; /* accumulate ops may leave values outside (-q,q) */
  if (c < 0)
    c += q;
  return c;
}

/* ------------------------------------------------------------------ */
/* vec                                                                */
/* ------------------------------------------------------------------ */

static arith_vec_ptr
lz_vec_alloc (const arith_ring_ptr r, unsigned int n)
{
  arith_vec_ptr v = (arith_vec_ptr)malloc (sizeof (struct arith_vec));
  v->ring = r;
  polyvec_alloc (v->v, r->ring, n);
  for (unsigned int i = 0; i < n; i++)
    poly_set_zero (polyvec_get_elem (v->v, i));
  return v;
}

static void
lz_vec_free (arith_vec_ptr v)
{
  polyvec_free (v->v);
  free (v);
}

static void
lz_vec_set (arith_vec_ptr dst, const arith_vec_ptr src)
{
  polyvec_set (dst->v, src->v);
}

static void
lz_vec_set_coeff (arith_vec_ptr v, unsigned int elem, unsigned int k,
                  long long c)
{
  int_ptr coeff = poly_get_coeff (polyvec_get_elem (v->v, elem), k);
  int_set_i64 (coeff, c);
}

static long long
lz_vec_get_coeff (const arith_vec_ptr v, unsigned int elem, unsigned int k)
{
  long long c = int_get_i64 (
      poly_get_coeff (polyvec_get_elem (v->v, elem), k));
  long long q = (long long)v->ring->q;
  c %= q; /* accumulate ops may leave values outside (-q,q) */
  if (c < 0)
    c += q;
  return c;
}

/* ------------------------------------------------------------------ */
/* mat                                                                */
/* ------------------------------------------------------------------ */

static arith_mat_ptr
lz_mat_alloc (const arith_ring_ptr r, unsigned int nrows, unsigned int ncols)
{
  arith_mat_ptr A = (arith_mat_ptr)malloc (sizeof (struct arith_mat));
  A->ring = r;
  polymat_alloc (A->A, r->ring, nrows, ncols);
  for (unsigned int i = 0; i < nrows; i++)
    for (unsigned int j = 0; j < ncols; j++)
      poly_set_zero (polymat_get_elem (A->A, i, j));
  return A;
}

static void
lz_mat_free (arith_mat_ptr A)
{
  polymat_free (A->A);
  free (A);
}

static void
lz_mat_set_coeff (arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k, long long c)
{
  int_ptr coeff = poly_get_coeff (polymat_get_elem (A->A, i, j), k);
  int_set_i64 (coeff, c);
}

static long long
lz_mat_get_coeff (const arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k)
{
  long long c = int_get_i64 (
      poly_get_coeff (polymat_get_elem (A->A, i, j), k));
  long long q = (long long)A->ring->q;
  c %= q; /* accumulate ops may leave values outside (-q,q) */
  if (c < 0)
    c += q;
  return c;
}

/* ------------------------------------------------------------------ */
/* algebra                                                            */
/* ------------------------------------------------------------------ */

static void
lz_vec_addmul (const arith_vec_ptr w, const arith_mat_ptr A,
               const arith_vec_ptr v)
{
  polyvec_addmul (w->v, A->A, v->v, 0);
}

static void
lz_vec_acc_polyscalar (const arith_vec_ptr w, const arith_poly_ptr s,
                       const arith_vec_ptr v)
{
  polyvec_addscale2 (w->v, s->p, v->v, 0);
}

static void
lz_vec_acc (const arith_vec_ptr w, const arith_vec_ptr a, int neg)
{
  if (neg)
    polyvec_sub (w->v, w->v, a->v, 0);
  else
    polyvec_add (w->v, w->v, a->v, 0);
}

/* ------------------------------------------------------------------ */

static const arith_ops_s ops = {
    "lazer",

    lz_ring_new,
    lz_ring_free,
    lz_ring_q,
    lz_ring_d,

    lz_poly_alloc,
    lz_poly_free,
    lz_poly_set_coeff,
    lz_poly_get_coeff,

    lz_vec_alloc,
    lz_vec_free,
    lz_vec_set,
    lz_vec_set_coeff,
    lz_vec_get_coeff,

    lz_mat_alloc,
    lz_mat_free,
    lz_mat_set_coeff,
    lz_mat_get_coeff,

    lz_vec_addmul,
    lz_vec_acc_polyscalar,
    lz_vec_acc,
};

arith_ops_ptr
arith_lazer_ops (void)
{
  return &ops;
}