/**
 * \file   backend_ntl.cpp
 * \brief  NTL backend for the pluggable arithmetic interface.
 *
 * Ring elements are NTL ZZ_pX polynomials over Z/qZ reduced modulo
 * x^d + 1 (negacyclic wrap).  q must be prime (NTL ZZ_p requirement) and
 * below 2^63; ring construction from (q,d) alone with the optional `spec`
 * ignored.
 *
 * NTL arithmetic happens under an explicitly restored ZZ_pContext.  Every
 * backend operation opens a ctx_guard scope that pushes the ring's context
 * (make it current) and pops it again on exit, so backend calls nest safely
 * and never leak context-stack entries into the caller.
 */

#include "arith.h"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <cstdlib>

struct arith_ring
{
  NTL::ZZ_pContext *ctx;
  unsigned long long q;
  unsigned int d;
};

struct arith_poly
{
  struct arith_ring *r;
  NTL::ZZ_pX p;
};

struct arith_vec
{
  struct arith_ring *r;
  unsigned int n;
  NTL::ZZ_pX *e; /* n elements */
};

struct arith_mat
{
  struct arith_ring *r;
  unsigned int nr, nc;
  NTL::ZZ_pX *e; /* nr * nc elements, row-major */
};

namespace
{
/* Make the ring's NTL context current for the lifetime of the scope. */
struct ctx_guard
{
  NTL::ZZ_pBak bak;
  explicit ctx_guard (const struct arith_ring *r)
  {
    bak.save ();
    r->ctx->restore ();
  }
  ~ctx_guard ()
  {
    bak.restore ();
  }
};

/* res = a * b mod (x^d + 1); NTL context must be current; res must not
 * alias a or b. */
static void
mul_negacyclic (const struct arith_ring *r, NTL::ZZ_pX &res,
                const NTL::ZZ_pX &a, const NTL::ZZ_pX &b)
{
  res = a * b;
  long degOut = NTL::deg (res);
  long d = (long)r->d;
  if (degOut >= d)
    {
      for (long k = d; k <= degOut; k++)
        {
          NTL::ZZ_p ck = NTL::coeff (res, k);
          if (NTL::IsZero (ck))
            continue;
          NTL::ZZ_p ckmd = NTL::coeff (res, k - d);
          NTL::SetCoeff (res, k - d, ckmd - ck);
        }
      for (long k = d; k <= degOut; k++)
        NTL::SetCoeff (res, k, 0);
      res.normalize ();
    }
}
}

/* ------------------------------------------------------------------ */
/* ring                                                                */
/* ------------------------------------------------------------------ */

static arith_ring_ptr
nt_ring_new (unsigned long long q, unsigned int d, const void *spec)
{
  (void)spec;

  if (q == 0 || q >= (1ULL << 63))
    return NULL;

  struct arith_ring *r = new struct arith_ring;
  r->ctx = new NTL::ZZ_pContext (NTL::to_ZZ (static_cast<long>(q)));
  r->q = q;
  r->d = d;
  return r;
}

static void
nt_ring_free (arith_ring_ptr r)
{
  delete r->ctx;
  delete r;
}

static unsigned long long
nt_ring_q (const arith_ring_ptr r)
{
  return r->q;
}

static unsigned int
nt_ring_d (const arith_ring_ptr r)
{
  return r->d;
}

/* ------------------------------------------------------------------ */
/* poly                                                                */
/* ------------------------------------------------------------------ */

static arith_poly_ptr
nt_poly_alloc (const arith_ring_ptr r)
{
  arith_poly_ptr p = new struct arith_poly;
  p->r = r;
  return p;
}

static void
nt_poly_free (arith_poly_ptr p)
{
  delete p;
}

static void
nt_poly_set_coeff (arith_poly_ptr p, unsigned int i, long long c)
{
  ctx_guard g (p->r);
  c %= (long long)p->r->q;
  if (c < 0)
    c += (long long)p->r->q;
  NTL::SetCoeff (p->p, (long)i, NTL::conv<NTL::ZZ_p> (c));
}

static long long
nt_poly_get_coeff (const arith_poly_ptr p, unsigned int i)
{
  ctx_guard g (p->r);
  return (long long)NTL::conv<long> (NTL::rep (NTL::coeff (p->p, (long)i)));
}

/* ------------------------------------------------------------------ */
/* vec                                                                */
/* ------------------------------------------------------------------ */

static arith_vec_ptr
nt_vec_alloc (const arith_ring_ptr r, unsigned int n)
{
  arith_vec_ptr v = new struct arith_vec;
  v->r = r;
  v->n = n;
  v->e = new NTL::ZZ_pX[n];
  return v;
}

static void
nt_vec_free (arith_vec_ptr v)
{
  delete[] v->e;
  delete v;
}

static void
nt_vec_set (arith_vec_ptr dst, const arith_vec_ptr src)
{
  for (unsigned int i = 0; i < dst->n; i++)
    dst->e[i] = src->e[i];
}

static void
nt_vec_set_coeff (arith_vec_ptr v, unsigned int elem, unsigned int k,
                  long long c)
{
  ctx_guard g (v->r);
  c %= (long long)v->r->q;
  if (c < 0)
    c += (long long)v->r->q;
  NTL::SetCoeff (v->e[elem], (long)k, NTL::conv<NTL::ZZ_p> (c));
}

static long long
nt_vec_get_coeff (const arith_vec_ptr v, unsigned int elem, unsigned int k)
{
  ctx_guard g (v->r);
  return (long long)NTL::conv<long> (
      NTL::rep (NTL::coeff (v->e[elem], (long)k)));
}

/* ------------------------------------------------------------------ */
/* mat                                                                */
/* ------------------------------------------------------------------ */

static arith_mat_ptr
nt_mat_alloc (const arith_ring_ptr r, unsigned int nrows, unsigned int ncols)
{
  arith_mat_ptr A = new struct arith_mat;
  A->r = r;
  A->nr = nrows;
  A->nc = ncols;
  A->e = new NTL::ZZ_pX[nrows * ncols];
  return A;
}

static void
nt_mat_free (arith_mat_ptr A)
{
  delete[] A->e;
  delete A;
}

static void
nt_mat_set_coeff (arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k, long long c)
{
  ctx_guard g (A->r);
  c %= (long long)A->r->q;
  if (c < 0)
    c += (long long)A->r->q;
  NTL::SetCoeff (A->e[i * A->nc + j], (long)k, NTL::conv<NTL::ZZ_p> (c));
}

static long long
nt_mat_get_coeff (const arith_mat_ptr A, unsigned int i, unsigned int j,
                  unsigned int k)
{
  ctx_guard g (A->r);
  return (long long)NTL::conv<long> (
      NTL::rep (NTL::coeff (A->e[i * A->nc + j], (long)k)));
}

/* ------------------------------------------------------------------ */
/* algebra                                                            */
/* ------------------------------------------------------------------ */

static void
nt_vec_addmul (const arith_vec_ptr w, const arith_mat_ptr A,
               const arith_vec_ptr v)
{
  ctx_guard g (w->r);

  for (unsigned int i = 0; i < A->nr; i++)
    {
      NTL::ZZ_pX acc;
      for (unsigned int j = 0; j < A->nc; j++)
        {
          NTL::ZZ_pX t;
          mul_negacyclic (w->r, t, A->e[i * A->nc + j], v->e[j]);
          acc += t;
        }
      w->e[i] += acc;
    }
}

static void
nt_vec_acc_polyscalar (const arith_vec_ptr w, const arith_poly_ptr s,
                       const arith_vec_ptr v)
{
  ctx_guard g (w->r);

  for (unsigned int i = 0; i < w->n; i++)
    {
      NTL::ZZ_pX t;
      mul_negacyclic (w->r, t, s->p, v->e[i]);
      w->e[i] += t;
    }
}

static void
nt_vec_acc (const arith_vec_ptr w, const arith_vec_ptr a, int neg)
{
  ctx_guard g (w->r);

  for (unsigned int i = 0; i < w->n; i++)
    {
      if (neg)
        w->e[i] -= a->e[i];
      else
        w->e[i] += a->e[i];
    }
}

/* ------------------------------------------------------------------ */

static const arith_ops_s ops = {
    "ntl",

    nt_ring_new,
    nt_ring_free,
    nt_ring_q,
    nt_ring_d,

    nt_poly_alloc,
    nt_poly_free,
    nt_poly_set_coeff,
    nt_poly_get_coeff,

    nt_vec_alloc,
    nt_vec_free,
    nt_vec_set,
    nt_vec_set_coeff,
    nt_vec_get_coeff,

    nt_mat_alloc,
    nt_mat_free,
    nt_mat_set_coeff,
    nt_mat_get_coeff,

    nt_vec_addmul,
    nt_vec_acc_polyscalar,
    nt_vec_acc,
};

extern "C" arith_ops_ptr
arith_ntl_ops (void)
{
  return &ops;
}