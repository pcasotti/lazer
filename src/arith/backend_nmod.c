/**
 * \file   backend_nmod.c
 * \brief  FLINT nmod_poly backend for the pluggable arithmetic interface.
 *
 * Ring elements are nmod_poly's over Z/qZ reduced modulo x^d + 1 (negacyclic
 * wrap).  Single-limb only: q must fit in one machine word, and (per the
 * interface's long long accessor contract) must stay below 2^63.  Ring
 * construction from (q,d) alone; the optional `spec` is ignored.
 */

#include "arith.h"

#include <flint/nmod.h>
#include <flint/nmod_poly.h>

#include <stdlib.h>

struct arith_ring
{
  nmod_t mod; /* mod.n == q */
  unsigned long long q;
  unsigned int d;
};

struct arith_poly
{
  struct arith_ring *r;
  nmod_poly_struct p;
};

struct arith_vec
{
  struct arith_ring *r;
  unsigned int n;
  nmod_poly_struct *e; /* n elements */
};

struct arith_mat
{
  struct arith_ring *r;
  unsigned int nr, nc;
  nmod_poly_struct *e; /* nr * nc elements, row-major */
};

static inline nmod_poly_struct *
mp_elem (const struct arith_mat *A, unsigned int i, unsigned int j)
{
  return (nmod_poly_struct *)&A->e[i * A->nc + j];
}

static inline nmod_t
mod_of (const struct arith_ring *r)
{
  return r->mod;
}

/* res = a * b mod (x^d + 1); res must not alias a or b */
static void
mul_negacyclic (nmod_poly_struct *res, const nmod_poly_struct *a,
                const nmod_poly_struct *b, unsigned int d)
{
  nmod_poly_mul (res, a, b);

  slong deg = nmod_poly_degree (res);
  if (deg >= (slong)d)
    {
      nmod_t mod = res->mod;
      for (slong k = d; k <= deg; k++)
        {
          ulong ck = nmod_poly_get_coeff_ui (res, k);
          if (ck == 0)
            continue;
          ulong cur = nmod_poly_get_coeff_ui (res, k - d);
          nmod_poly_set_coeff_ui (res, k - d, nmod_sub (cur, ck, mod));
        }
      nmod_poly_truncate (res, d);
    }
}

/* ------------------------------------------------------------------ */
/* ring                                                                */
/* ------------------------------------------------------------------ */

static arith_ring_ptr
np_ring_new (unsigned long long q, unsigned int d, const void *spec)
{
  (void)spec;

  if (q == 0 || q >= (1ULL << 63))
    return NULL; /* single limb, and [0,q) must fit long long */

  struct arith_ring *r = (struct arith_ring *)malloc (sizeof (struct arith_ring));
  nmod_init (&r->mod, q);
  r->q = q;
  r->d = d;
  return r;
}

static void
np_ring_free (arith_ring_ptr r)
{
  free (r);
}

static unsigned long long
np_ring_q (const arith_ring_ptr r)
{
  return r->q;
}

static unsigned int
np_ring_d (const arith_ring_ptr r)
{
  return r->d;
}

/* ------------------------------------------------------------------ */
/* poly                                                                */
/* ------------------------------------------------------------------ */

static arith_poly_ptr
np_poly_alloc (const arith_ring_ptr r)
{
  arith_poly_ptr p = (arith_poly_ptr)malloc (sizeof (struct arith_poly));
  p->r = r;
  nmod_poly_init (&p->p, r->q);
  return p;
}

static void
np_poly_free (arith_poly_ptr p)
{
  nmod_poly_clear (&p->p);
  free (p);
}

static void
np_poly_set_coeff (arith_poly_ptr p, unsigned int i, long long c)
{
  c %= (long long)mod_of (p->r).n;
  if (c < 0)
    c += (long long)mod_of (p->r).n;
  nmod_poly_set_coeff_ui (&p->p, i, (ulong)c);
}

static long long
np_poly_get_coeff (const arith_poly_ptr p, unsigned int i)
{
  return (long long)nmod_poly_get_coeff_ui (&p->p, i);
}

/* ------------------------------------------------------------------ */
/* vec                                                                */
/* ------------------------------------------------------------------ */

static arith_vec_ptr
np_vec_alloc (const arith_ring_ptr r, unsigned int n)
{
  arith_vec_ptr v = (arith_vec_ptr)malloc (sizeof (struct arith_vec));
  v->r = r;
  v->n = n;
  v->e = flint_malloc (n * sizeof (nmod_poly_struct));
  for (unsigned int i = 0; i < n; i++)
    nmod_poly_init (&v->e[i], r->q);
  return v;
}

static void
np_vec_free (arith_vec_ptr v)
{
  for (unsigned int i = 0; i < v->n; i++)
    nmod_poly_clear (&v->e[i]);
  flint_free (v->e);
  free (v);
}

static void
np_vec_set (arith_vec_ptr dst, const arith_vec_ptr src)
{
  for (unsigned int i = 0; i < dst->n; i++)
    nmod_poly_set (&dst->e[i], &src->e[i]);
}

static void
np_vec_set_coeff (arith_vec_ptr v, unsigned int elem, unsigned int k,
                  long long c)
{
  c %= (long long)mod_of (v->r).n;
  if (c < 0)
    c += (long long)mod_of (v->r).n;
  nmod_poly_set_coeff_ui (&v->e[elem], k, (ulong)c);
}

static long long
np_vec_get_coeff (const arith_vec_ptr v, unsigned int elem, unsigned int k)
{
  return (long long)nmod_poly_get_coeff_ui (&v->e[elem], k);
}

/* ------------------------------------------------------------------ */
/* mat                                                                */
/* ------------------------------------------------------------------ */

static arith_mat_ptr
np_mat_alloc (const arith_ring_ptr r, unsigned int nrows, unsigned int ncols)
{
  arith_mat_ptr A = (arith_mat_ptr)malloc (sizeof (struct arith_mat));
  A->r = r;
  A->nr = nrows;
  A->nc = ncols;
  A->e = flint_malloc (nrows * ncols * sizeof (nmod_poly_struct));
  for (unsigned int i = 0; i < nrows * ncols; i++)
    nmod_poly_init (&A->e[i], r->q);
  return A;
}

static void
np_mat_free (arith_mat_ptr A)
{
  for (unsigned int i = 0; i < A->nr * A->nc; i++)
    nmod_poly_clear (&A->e[i]);
  flint_free (A->e);
  free (A);
}

static void
np_mat_set_coeff (arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k, long long c)
{
  c %= (long long)mod_of (A->r).n;
  if (c < 0)
    c += (long long)mod_of (A->r).n;
  nmod_poly_set_coeff_ui (mp_elem (A, i, j), k, (ulong)c);
}

static long long
np_mat_get_coeff (const arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k)
{
  return (long long)nmod_poly_get_coeff_ui (mp_elem (A, i, j), k);
}

/* ------------------------------------------------------------------ */
/* algebra                                                            */
/* ------------------------------------------------------------------ */

static void
np_vec_addmul (const arith_vec_ptr w, const arith_mat_ptr A,
               const arith_vec_ptr v)
{
  unsigned int d = w->r->d;

  nmod_poly_struct *acc = flint_malloc (sizeof (nmod_poly_struct));
  nmod_poly_struct *t = flint_malloc (sizeof (nmod_poly_struct));
  nmod_poly_init (acc, w->r->q);
  nmod_poly_init (t, w->r->q);

  for (unsigned int i = 0; i < A->nr; i++)
    {
      nmod_poly_zero (acc);
      for (unsigned int j = 0; j < A->nc; j++)
        {
          mul_negacyclic (t, mp_elem (A, i, j), &v->e[j], d);
          nmod_poly_add (acc, acc, t);
        }
      nmod_poly_add (&w->e[i], &w->e[i], acc);
    }

  nmod_poly_clear (t);
  nmod_poly_clear (acc);
  flint_free (t);
  flint_free (acc);
}

static void
np_vec_acc_polyscalar (const arith_vec_ptr w, const arith_poly_ptr s,
                       const arith_vec_ptr v)
{
  unsigned int d = w->r->d;

  nmod_poly_struct *t = flint_malloc (sizeof (nmod_poly_struct));
  nmod_poly_init (t, w->r->q);

  for (unsigned int i = 0; i < w->n; i++)
    {
      mul_negacyclic (t, &s->p, &v->e[i], d);
      nmod_poly_add (&w->e[i], &w->e[i], t);
    }

  nmod_poly_clear (t);
  flint_free (t);
}

static void
np_vec_acc (const arith_vec_ptr w, const arith_vec_ptr a, int neg)
{
  for (unsigned int i = 0; i < w->n; i++)
    {
      if (neg)
        nmod_poly_sub (&w->e[i], &w->e[i], &a->e[i]);
      else
        nmod_poly_add (&w->e[i], &w->e[i], &a->e[i]);
    }
}

/* ------------------------------------------------------------------ */

static const arith_ops_s ops = {
    "nmod",

    np_ring_new,
    np_ring_free,
    np_ring_q,
    np_ring_d,

    np_poly_alloc,
    np_poly_free,
    np_poly_set_coeff,
    np_poly_get_coeff,

    np_vec_alloc,
    np_vec_free,
    np_vec_set,
    np_vec_set_coeff,
    np_vec_get_coeff,

    np_mat_alloc,
    np_mat_free,
    np_mat_set_coeff,
    np_mat_get_coeff,

    np_vec_addmul,
    np_vec_acc_polyscalar,
    np_vec_acc,
};

arith_ops_ptr
arith_nmod_ops (void)
{
  return &ops;
}