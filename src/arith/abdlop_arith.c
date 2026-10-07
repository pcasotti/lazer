/*
 * \file   abdlop_arith.c
 * \brief  abdlop pipeline over the pluggable arith interface.
 *
 * Mirrors src/abdlop.c (keygen, commit, hashcomm, prove, verify) with every
 * ring multiplication going through the arith ops; compression, encoding,
 * hashing, gaussian/autostable sampling and rejection stay above the
 * interface, coefficient by coefficient.  Deterministic input generation
 * (uniform/short samplers on the seed) also lives above the interface.
 */

#include "abdlop_arith.h"

#include <stdint.h>
#include <string.h>
#include <time.h>

/* ========================= coeff helpers ========================= */

static long long
_canonq (long long v, unsigned long long q)
{
  if (v < 0)
    v += (long long)q;
  else if (v >= (long long)q)
    v -= (long long)q;
  return v;
}

static long long
_balq (long long c, unsigned long long q)
{
  const long long pivot = (long long)((q - 1) >> 1);
  return c > pivot ? c - (long long)q : c;
}

/* decomp_c: canonical c in [0,q) -> (r1, r0), mirrors dcompress_decompose. */
static void
_decomp_c (long long c, long long *r1p, long long *r0p,
           const abdlop_arith_params_t *ap)
{
  const long long gamma = ap->gamma, gb2 = ap->gammaby2;
  long long r0 = c % gamma;
  long long r1;

  if (r0 > gb2)
    r0 -= gamma;
  r1 = c - r0;
  if (r1 == ap->qminus1)
    {
      r1 = 0;
      r0 -= 1;
    }
  else
    r1 /= gamma;
  *r1p = r1;
  *r0p = r0;
}

/* p2r_c: canonical c in [0,q) -> (t1,t2) = power2round, mirrors dcompress. */
static long long
_p2r_c (long long c, long long *t2p, const abdlop_arith_params_t *ap)
{
  const long long pd = ap->pow2D, pd2 = ap->pow2Dby2;
  long long r0 = c & (pd - 1);
  long long t1, t2;

  if (r0 > pd2)
    r0 -= pd;
  t1 = (c - r0) >> ap->D;
  t2 = c - (t1 << ap->D);
  *t2p = _canonq (t2, ap->q);
  return t1;
}

/* use_ghint_c: r canonical in [0,q), y signed hint -> result in [0,m). */
static long long
_use_ghint_c (long long r, long long y, const abdlop_arith_params_t *ap)
{
  long long r1, r0, s;

  _decomp_c (r, &r1, &r0, ap);
  s = r1 + y;
  s %= ap->m;
  if (s < 0)
    s += ap->m;
  return s;
}

/* make_ghint_c: r canonical in [0,q), z signed -> hint, mirrors dcompress. */
static long long
_make_ghint_c (long long r, long long z, const abdlop_arith_params_t *ap)
{
  const long long m = ap->m, mb2 = ap->mby2;
  long long s, r1a, r0a, r1b, r0b, h;

  _decomp_c (r, &r1a, &r0a, ap);
  s = _canonq (r + z, ap->q);
  _decomp_c (s, &r1b, &r0b, ap);
  h = r1b - r1a;
  h %= m;
  if (h < 0)
    h += m;
  if (h < 0)
    {
      if (-h >= mb2)
        h += m;
    }
  else
    {
      if (h > mb2)
        h -= m;
    }
  return h;
}

static __int128
_l2sqr (abdlop_arith_engine_t *e, const arith_vec_ptr v, unsigned int n,
        const abdlop_arith_params_t *ap)
{
  __int128 s = 0;
  unsigned int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      {
        long long b = _balq (e->ops->vec_get_coeff (v, i, k), ap->q);
        s += (__int128)b * b;
      }
  return s;
}

static long long
_linf (abdlop_arith_engine_t *e, const arith_vec_ptr v, unsigned int n,
       const abdlop_arith_params_t *ap)
{
  long long mx = 0;
  unsigned int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      {
        long long b = _balq (e->ops->vec_get_coeff (v, i, k), ap->q);
        if (b < 0)
          b = -b;
        if (b > mx)
          mx = b;
      }
  return mx;
}

/* ==================== above-interface mirrors ==================== */

/* copy canonical coefficients of an interface vec into a lazer polyvec. */
static void
_v2lz (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
       polyvec_t lz, const arith_vec_ptr iv, unsigned int n)
{
  unsigned int i, k;
  INT_T (sc, int_get_nlimbs (polyring_get_mod (polyvec_get_ring (lz))));

  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      {
        int_set_i64 (sc, e->ops->vec_get_coeff (iv, i, k));
        poly_set_coeff (polyvec_get_elem (lz, i), k, sc);
      }
}

/* copy signed (balanced) coefficients of an interface vec into an intvec. */
static void
_v2int (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
        intvec_t r, const arith_vec_ptr iv, unsigned int n)
{
  unsigned int i, k;
  unsigned int idx = 0;

  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      intvec_set_elem_i64 (r, idx++, _balq (e->ops->vec_get_coeff (iv, i, k),
                                            ap->q));
}

/* import a lazer polyvec into an interface vec through the accessors. */
static void
_lz2v (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
       arith_vec_ptr iv, polyvec_t lz, unsigned int n)
{
  unsigned int i, k;

  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      e->ops->vec_set_coeff (
          iv, i, k,
          _canonq (int_get_i64 (poly_get_coeff (polyvec_get_elem (lz, i), k)),
                   ap->q));
}

/* encode an interface vec like coder_enc_urandom3 over a fresh mirror. */
static void
_enc_urandom (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
              coder_state_t cstate, const arith_vec_ptr iv, unsigned int n,
              const int_t mod, unsigned int mbits)
{
  polyvec_t mrc;

  polyvec_alloc (mrc, ap->p->ring, n);
  _v2lz (e, ap, mrc, iv, n);
  polyvec_redp (mrc, mrc);
  coder_enc_urandom3 (cstate, mrc, mod, mbits);
  polyvec_free (mrc);
}

/* ========================= pipeline ========================= */

void
abdlop_arith_params_init (abdlop_arith_params_t *ap,
                          const abdlop_params_srcptr params)
{
  const dcompress_params_srcptr d = params->dcompress;

  ap->q = int_get_i64 (params->ring->q);
  ap->d = polyring_get_deg (params->ring);
  ap->log2q = polyring_get_log2q (params->ring);
  ap->D = d->D;
  ap->log2m = d->log2m;
  ap->log2qmd = ap->log2q - ap->D;
  ap->m1 = params->m1;
  ap->m2 = params->m2;
  ap->m22 = params->m2 - params->kmsis;
  ap->kmsis = params->kmsis;
  ap->l = params->l + params->lext;
  ap->qminus1 = int_get_i64 (d->qminus1);
  ap->m = int_get_i64 (d->m);
  ap->mby2 = int_get_i64 (d->mby2);
  ap->gamma = int_get_i64 (d->gamma);
  ap->gammaby2 = int_get_i64 (d->gammaby2);
  ap->pow2D = int_get_i64 (d->pow2D);
  ap->pow2Dby2 = int_get_i64 (d->pow2Dby2);
  ap->Bsqr = int_get_i64 (params->Bsqr);
  ap->m_odd = d->m_odd;
  ap->omega = (unsigned int)params->omega;
  ap->log2omega = params->log2omega;
  ap->log2stdev1 = params->log2stdev1;
  ap->log2stdev2 = params->log2stdev2;
  ap->rej1 = params->rej1;
  ap->rej2 = params->rej2;
  ap->p = params;
}

void
abdlop_arith_keygen (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
                     const uint8_t seed[32], arith_mat_ptr A1,
                     arith_mat_ptr A2prime, arith_mat_ptr Bprime)
{
  polyring_srcptr Rq = ap->p->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int m1 = ap->m1;
  const unsigned int m22 = ap->m22;
  const unsigned int l = ap->l;
  const unsigned int kmsis = ap->kmsis;
  polymat_t a1, a2, b;
  unsigned int i, j, k;

  if (m1 > 0)
    {
      polymat_alloc (a1, Rq, kmsis, m1);
      polymat_alloc (a2, Rq, kmsis, m22);
      polymat_urandom (a1, q, log2q, seed, 0);
      polymat_urandom (a2, q, log2q, seed, 1);
      for (i = 0; i < kmsis; i++)
        for (j = 0; j < m1; j++)
          for (k = 0; k < ap->d; k++)
            e->ops->mat_set_coeff (
                A1, i, j, k,
                _canonq (int_get_i64 (poly_get_coeff (
                            polymat_get_elem (a1, i, j), k)),
                         ap->q));
      for (i = 0; i < kmsis; i++)
        for (j = 0; j < m22; j++)
          for (k = 0; k < ap->d; k++)
            e->ops->mat_set_coeff (
                A2prime, i, j, k,
                _canonq (int_get_i64 (poly_get_coeff (
                            polymat_get_elem (a2, i, j), k)),
                         ap->q));
      polymat_free (a1);
      polymat_free (a2);
    }

  if (l > 0)
    {
      polymat_alloc (b, Rq, l, m22);
      polymat_urandom (b, q, log2q, seed, 2);
      for (i = 0; i < l; i++)
        for (j = 0; j < m22; j++)
          for (k = 0; k < ap->d; k++)
            e->ops->mat_set_coeff (
                Bprime, i, j, k,
                _canonq (int_get_i64 (poly_get_coeff (
                            polymat_get_elem (b, i, j), k)),
                         ap->q));
      polymat_free (b);
    }
}

void
abdlop_arith_commit (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
                     arith_vec_ptr tA1, arith_vec_ptr tA2, arith_vec_ptr tB,
                     const arith_vec_ptr s1, const arith_vec_ptr m,
                     const arith_vec_ptr s2, const arith_mat_ptr A1,
                     const arith_mat_ptr A2prime, const arith_mat_ptr Bprime)
{
  const unsigned int m1 = ap->m1;
  const unsigned int m22 = ap->m22;
  const unsigned int kmsis = ap->kmsis;
  const unsigned int l = ap->l;
  const unsigned int d = ap->d;
  unsigned int i, k;
  arith_vec_ptr s21;

  /* split s2 = (s21, s22) */
  s21 = e->ops->vec_alloc (e->ring, m22);
  for (i = 0; i < m22; i++)
    for (k = 0; k < d; k++)
      e->ops->vec_set_coeff (s21, i, k,
                             e->ops->vec_get_coeff (s2, i, k));

  if (m1 > 0)
    {
      arith_vec_ptr s22 = e->ops->vec_alloc (e->ring, kmsis);
      for (i = 0; i < kmsis; i++)
        for (k = 0; k < d; k++)
          e->ops->vec_set_coeff (s22, i, k,
                                 e->ops->vec_get_coeff (s2, m22 + i, k));

      e->ops->vec_set (tA2, s22);
      e->ops->vec_addmul (tA2, A1, s1);
      e->ops->vec_addmul (tA2, A2prime, s21);

      for (i = 0; i < kmsis; i++)
        for (k = 0; k < d; k++)
          {
            long long c = e->ops->vec_get_coeff (tA2, i, k);
            long long t1v, t2v;
            t1v = _p2r_c (c, &t2v, ap);
            e->ops->vec_set_coeff (tA1, i, k, t1v);
            e->ops->vec_set_coeff (tA2, i, k, t2v);
          }
      e->ops->vec_free (s22);
    }

  if (l > 0)
    {
      e->ops->vec_set (tB, m);
      e->ops->vec_addmul (tB, Bprime, s21);
    }

  e->ops->vec_free (s21);
}

void
abdlop_arith_hashcomm (uint8_t hash[32], abdlop_arith_engine_t *e,
                       const abdlop_arith_params_t *ap,
                       const arith_vec_ptr tA1, const arith_vec_ptr tB)
{
  const unsigned int d = ap->d;
  const unsigned int kmsis = ap->kmsis;
  const unsigned int l = ap->l;
  coder_state_t cstate;
  shake128_state_t hstate;
  const size_t outlen
      = CEIL (kmsis * d * ap->log2qmd + l * d * ap->log2q, 8) + 1;
  uint8_t out[outlen];
  unsigned int nbits;

  coder_enc_begin (cstate, out);
  if (ap->m1 > 0)
    {
      INT_T (mod, int_get_nlimbs (polyring_get_mod (ap->p->ring)));
      int_set_one (mod);
      int_lshift (mod, mod, ap->log2qmd);
      _enc_urandom (e, ap, cstate, tA1, kmsis, mod, ap->log2qmd);
    }
  if (l > 0)
    _enc_urandom (e, ap, cstate, tB, l, ap->p->ring->q, ap->log2q);
  coder_enc_end (cstate);
  nbits = coder_get_offset (cstate) >> 3;

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, nbits);
  shake128_squeeze (hstate, hash, 32);
  shake128_clear (hstate);
}

/* challenge <- poly_urandom_autostable(H(hash || encode(w1))). */
static void
_challenge (uint8_t hash[32], abdlop_arith_engine_t *e,
            const abdlop_arith_params_t *ap, const arith_vec_ptr w1,
            poly_t clz, uint8_t cseed[32])
{
  const unsigned int d = ap->d;
  const unsigned int kmsis = ap->kmsis;
  coder_state_t cstate;
  shake128_state_t hstate;
  uint8_t out[CEIL (ap->log2m * d * kmsis, 8) + 1];
  uint8_t cseed_tmp[32];
  unsigned int nbits;

  coder_enc_begin (cstate, out);
  _enc_urandom (e, ap, cstate, w1, kmsis, ap->p->dcompress->m, ap->log2m);
  coder_enc_end (cstate);
  nbits = coder_get_offset (cstate) >> 3;

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, nbits);
  shake128_squeeze (hstate, cseed_tmp, 32);
  shake128_clear (hstate);

  poly_urandom_autostable (clz, (int64_t)ap->omega, ap->log2omega,
                           cseed_tmp, 0);
  memcpy (cseed, cseed_tmp, 32);
}

void
abdlop_arith_prove (uint8_t hash[32], arith_poly_ptr c, arith_vec_ptr z1,
                    arith_vec_ptr z21, arith_vec_ptr h,
                    const arith_vec_ptr tA2, const arith_vec_ptr s1,
                    const arith_vec_ptr s2, const arith_mat_ptr A1,
                    const arith_mat_ptr A2prime, const uint8_t seed[32],
                    abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
                    abdlop_arith_prove_stats_t *stats)
{
  const unsigned int m1 = ap->m1;
  const unsigned int m2 = ap->m2;
  const unsigned int m22 = ap->m22;
  const unsigned int kmsis = ap->kmsis;
  const unsigned int d = ap->d;
  const unsigned int nlimbs = int_get_nlimbs (polyring_get_mod (ap->p->ring));
  polyring_srcptr Rq = ap->p->ring;
  rng_state_t rngstate;
  uint8_t yseed[32];
  uint32_t dom;
  uint8_t cseed[32];
  arith_poly_ptr negc;
  arith_vec_ptr s21, s22;
  arith_vec_ptr y1raw, y21raw, y22raw, cs1, cs2, zerom1, zerom2;
  arith_vec_ptr w, w1, w0, w1post;
  intvec_t z1c, cs1c, z2c, cs2c;
  polyvec_t y1_m, y2_m;
  poly_t c_m;
  unsigned int i, k, attempts = 0;
  int rej;
  clock_t t_loop, t_post;

  /* split s2 */
  s21 = e->ops->vec_alloc (e->ring, m22);
  s22 = e->ops->vec_alloc (e->ring, kmsis);
  for (i = 0; i < m22; i++)
    for (k = 0; k < d; k++)
      e->ops->vec_set_coeff (s21, i, k,
                             e->ops->vec_get_coeff (s2, i, k));
  for (i = 0; i < kmsis; i++)
    for (k = 0; k < d; k++)
      e->ops->vec_set_coeff (s22, i, k,
                             e->ops->vec_get_coeff (s2, m22 + i, k));

  negc = e->ops->poly_alloc (e->ring);
  y1raw = e->ops->vec_alloc (e->ring, m1);
  y21raw = e->ops->vec_alloc (e->ring, m22);
  y22raw = e->ops->vec_alloc (e->ring, kmsis);
  cs1 = e->ops->vec_alloc (e->ring, m1);
  cs2 = e->ops->vec_alloc (e->ring, m2);
  zerom1 = e->ops->vec_alloc (e->ring, m1);
  zerom2 = e->ops->vec_alloc (e->ring, m2);
  w = e->ops->vec_alloc (e->ring, kmsis);
  w1 = e->ops->vec_alloc (e->ring, kmsis);
  w0 = e->ops->vec_alloc (e->ring, kmsis);
  w1post = e->ops->vec_alloc (e->ring, kmsis);

  intvec_alloc (z1c, d * m1, nlimbs);
  intvec_alloc (cs1c, d * m1, nlimbs);
  intvec_alloc (z2c, d * m2, nlimbs);
  intvec_alloc (cs2c, d * m2, nlimbs);

  polyvec_alloc (y1_m, Rq, m1);
  polyvec_alloc (y2_m, Rq, m2);
  poly_alloc (c_m, Rq);

  rng_init (rngstate, seed, 0);
  rng_urandom (rngstate, yseed, 32);

  dom = 0;
  t_loop = stats != NULL ? clock () : 0;
  while (1)
    {
      attempts++;
      polyvec_grandom (y1_m, ap->log2stdev1, yseed, dom);
      dom++;
      polyvec_grandom (y2_m, ap->log2stdev2, yseed, dom);
      dom++;

      _lz2v (e, ap, y1raw, y1_m, m1);
      {
        polyvec_t y21_m, y22_m;
        polyvec_get_subvec (y21_m, y2_m, 0, m22, 1);
        polyvec_get_subvec (y22_m, y2_m, m22, kmsis, 1);
        _lz2v (e, ap, y21raw, y21_m, m22);
        _lz2v (e, ap, y22raw, y22_m, kmsis);
      }

      /* w = y22 + A1*y1 + A2prime*y21 */
      e->ops->vec_set (w, y22raw);
      e->ops->vec_addmul (w, A1, y1raw);
      e->ops->vec_addmul (w, A2prime, y21raw);

      /* w1, w0 <- decomp(w) */
      for (i = 0; i < kmsis; i++)
        for (k = 0; k < d; k++)
          {
            long long r1, r0;
            _decomp_c (e->ops->vec_get_coeff (w, i, k), &r1, &r0, ap);
            e->ops->vec_set_coeff (w1, i, k, r1);
            e->ops->vec_set_coeff (w0, i, k, r0);
          }

      /* challenge from w1 */
      _challenge (hash, e, ap, w1, c_m, cseed);
      for (k = 0; k < d; k++)
        e->ops->poly_set_coeff (
            c, k, _canonq (int_get_i64 (poly_get_coeff (c_m, k)), ap->q));

      /* z1 = y1 + c*s1 ; y21 <- y21 + c*s21 ; y22 <- y22 + c*s22 */
      e->ops->vec_set (z1, y1raw);
      e->ops->vec_acc_polyscalar (z1, c, s1);
      e->ops->vec_acc_polyscalar (y21raw, c, s21);
      e->ops->vec_acc_polyscalar (y22raw, c, s22);

      /* cs1 = c*s1 ; cs2 = c*s2 (on zero vectors) */
      e->ops->vec_set (cs1, zerom1);
      e->ops->vec_acc_polyscalar (cs1, c, s1);
      e->ops->vec_set (cs2, zerom2);
      e->ops->vec_acc_polyscalar (cs2, c, s2);

      /* rejection sampling (above the interface) */
      if (ap->rej1)
        {
          _v2int (e, ap, z1c, z1, m1);
          _v2int (e, ap, cs1c, cs1, m1);
          rej = rej_standard (rngstate, z1c, cs1c, ap->p->scM1,
                              ap->p->stdev1sqr);
          if (rej)
            continue;
        }
      if (ap->rej2)
        {
          _v2int (e, ap, z2c, y21raw, m22);
          for (i = 0; i < kmsis; i++)
            for (k = 0; k < d; k++)
              intvec_set_elem_i64 (
                  z2c, m22 * d + i * d + k,
                  _balq (e->ops->vec_get_coeff (y22raw, i, k), ap->q));
          _v2int (e, ap, cs2c, cs2, m2);
          rej = rej_standard (rngstate, z2c, cs2c, ap->p->scM2,
                              ap->p->stdev2sqr);
          if (rej)
            continue;
        }

      /* y22 <- y22 - c*tA2 - w0 */
      for (k = 0; k < d; k++)
        {
          long long ck = _balq (e->ops->poly_get_coeff (c, k), ap->q);
          e->ops->poly_set_coeff (negc, k, _canonq (-ck, ap->q));
        }
      e->ops->vec_acc_polyscalar (y22raw, negc, tA2);
      e->ops->vec_acc (y22raw, w0, -1);

      /* norm check on the full y2 = (y21, y22) */
      {
        __int128 nrm = _l2sqr (e, y21raw, m22, ap)
                       + _l2sqr (e, y22raw, kmsis, ap);
        if (nrm > (__int128)ap->Bsqr)
          continue;
      }

      break;
    }
  if (stats != NULL)
    stats->rej_ms = (double)(clock () - t_loop) / CLOCKS_PER_SEC;

  t_post = stats != NULL ? clock () : 0;

  /* w1post = gamma*w1 - y22 */
  {
    arith_poly_ptr gammapoly = e->ops->poly_alloc (e->ring);
    e->ops->poly_set_coeff (gammapoly, 0, ap->gamma);
    e->ops->vec_acc_polyscalar (w1post, gammapoly, w1);
    e->ops->vec_acc (w1post, y22raw, -1);
    e->ops->poly_free (gammapoly);
  }

  /* h = make_ghint(z = y22, r = w1post) */
  for (i = 0; i < kmsis; i++)
    for (k = 0; k < d; k++)
      {
        long long zv = _balq (e->ops->vec_get_coeff (y22raw, i, k), ap->q);
        long long rv = e->ops->vec_get_coeff (w1post, i, k);
        e->ops->vec_set_coeff (h, i, k, _make_ghint_c (rv, zv, ap));
      }

  /* z21 <- y21 (interface) */
  for (i = 0; i < m22; i++)
    for (k = 0; k < d; k++)
      e->ops->vec_set_coeff (z21, i, k,
                             e->ops->vec_get_coeff (y21raw, i, k));

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);

  /* cleanup */
  rng_clear (rngstate);
  e->ops->poly_free (negc);
  e->ops->vec_free (s21);
  e->ops->vec_free (s22);
  e->ops->vec_free (y1raw);
  e->ops->vec_free (y21raw);
  e->ops->vec_free (y22raw);
  e->ops->vec_free (cs1);
  e->ops->vec_free (cs2);
  e->ops->vec_free (zerom1);
  e->ops->vec_free (zerom2);
  e->ops->vec_free (w);
  e->ops->vec_free (w1);
  e->ops->vec_free (w0);
  e->ops->vec_free (w1post);
  intvec_free (z1c);
  intvec_free (cs1c);
  intvec_free (z2c);
  intvec_free (cs2c);
  polyvec_free (y1_m);
  polyvec_free (y2_m);
  poly_free (c_m);

  if (stats != NULL)
    {
      stats->attempts = attempts;
      stats->post_ms = (double)(clock () - t_post) / CLOCKS_PER_SEC;
    }
}

int
abdlop_arith_verify (uint8_t hash[32], const arith_poly_ptr c,
                     const arith_vec_ptr z1, const arith_vec_ptr z21,
                     const arith_vec_ptr h, const arith_vec_ptr tA1,
                     const arith_mat_ptr A1, const arith_mat_ptr A2prime,
                     abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap)
{
  const unsigned int m1 = ap->m1;
  const unsigned int kmsis = ap->kmsis;
  const unsigned int d = ap->d;
  arith_poly_ptr negshift;
  arith_vec_ptr tmp1, rem, w1v;
  poly_t c2;
  uint8_t cseed[32];
  unsigned int i, k;
  int b, accept = 0;

  negshift = e->ops->poly_alloc (e->ring);
  tmp1 = e->ops->vec_alloc (e->ring, kmsis);
  rem = e->ops->vec_alloc (e->ring, kmsis);
  w1v = e->ops->vec_alloc (e->ring, kmsis);

  /* -(c << D) coefficient-wise, mod q */
  for (k = 0; k < d; k++)
    {
      long long ck = _balq (e->ops->poly_get_coeff (c, k), ap->q);
      long long v = -(ck << ap->D);
      e->ops->poly_set_coeff (
          negshift, k,
          (v % (long long)ap->q + (long long)ap->q) % (long long)ap->q);
    }

  /* tmp1 = A1*z1 + A2prime*z21 - (c<<D)*tA1 */
  e->ops->vec_addmul (tmp1, A1, z1);
  e->ops->vec_addmul (tmp1, A2prime, z21);
  e->ops->vec_acc_polyscalar (tmp1, negshift, tA1);

  /* w1 <- use_ghint(h, tmp1) */
  for (i = 0; i < kmsis; i++)
    for (k = 0; k < d; k++)
      {
        long long yv = _balq (e->ops->vec_get_coeff (h, i, k), ap->q);
        long long rv = e->ops->vec_get_coeff (tmp1, i, k);
        e->ops->vec_set_coeff (w1v, i, k, _use_ghint_c (rv, yv, ap));
      }

  /* recover challenge from w1, compare to c */
  poly_alloc (c2, ap->p->ring);
  _challenge (hash, e, ap, w1v, c2, cseed);
  b = 1;
  for (k = 0; k < d; k++)
    {
      long long a = e->ops->poly_get_coeff (c, k);
      long long b2 = _canonq (int_get_i64 (poly_get_coeff (c2, k)), ap->q);
      if (a != b2)
        {
          b = 0;
          break;
        }
    }
  if (!b)
    goto ret;

  /* l2(z1)^2 <= sigma1^2 * 2 * m1 * d  (2.4025*400 = 962) */
  {
    __int128 bnd = (__int128)(2 * m1 * d * 962);
    bnd <<= 2 * ap->log2stdev1;
    bnd /= 400;
    b = _l2sqr (e, z1, m1, ap) <= bnd;
    if (!b)
      goto ret;
  }

  /* ||tmp1 - gamma*w1||^2 <= Bsqr */
  {
    arith_poly_ptr neggamma = e->ops->poly_alloc (e->ring);
    e->ops->poly_set_coeff (neggamma, 0, _canonq (-ap->gamma, ap->q));
    e->ops->vec_set (rem, tmp1);
    e->ops->vec_acc_polyscalar (rem, neggamma, w1v);
    b = _l2sqr (e, rem, kmsis, ap) <= (__int128)ap->Bsqr;
    e->ops->poly_free (neggamma);
    if (!b)
      goto ret;
  }

  /* linf(h) <= m */
  b = _linf (e, h, kmsis, ap) <= ap->m;
  if (!b)
    goto ret;

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);
  accept = 1;
ret:
  poly_free (c2);
  e->ops->poly_free (negshift);
  e->ops->vec_free (tmp1);
  e->ops->vec_free (rem);
  e->ops->vec_free (w1v);
  return accept;
}
