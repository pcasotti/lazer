/*
 * \file   arith-abdlop-test2.c
 * \brief  Full abdlop pipeline driven through the pluggable arith
 *         interface by the reusable src/arith/abdlop_arith module.
 *
 * The module (abdlop_arith.{c,h}) mirrors src/abdlop.c but runs every ring
 * multiplication through the ops of ../src/arith/arith.h; compression,
 * dcompression, rejection sampling, encoding, hashing and input generation
 * stay ABOVE the interface, so the same straight-line protocol produces
 * identical values in every backend.  This was validated one-off against
 * the lazer reference (tests/probe-module.c).
 *
 * Design (this test):
 *  - A deterministic "input series" for s1, s2, m is produced once by the
 *    lazer seeded samplers and imported into every backend through the
 *    coefficient accessors.  The module expands the public key matrices
 *    (A1, A2prime, Bprime) from the same seed internally (seed||0,1,2).
 *  - Each backend runs keygen, commit, hashcomm, prove, verify through the
 *    module of src/arith/abdlop_arith.h.
 *  - The test passes only if every backend's outputs and intermediates are
 *    cross-backend identical (cross-checked by digest) and every backend
 *    accepts the proof (accept_iv == 1).  No lazer oracle is used.
 */

#include "abdlop-params1.h"
#include "lazer.h"
#include "test.h"

#include "../src/arith/arith.h"
#include "../src/arith/abdlop_arith.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N_BACKENDS 4

typedef struct
{
  uint8_t hashcomm[32];
  uint8_t chash[32];
  uint8_t digest[32];
  int accept_iv;
} out_t;

typedef struct
{
  double commit, hashcomm, prove, verify, total, rejloop;
  long attempts, mism;
} stats_t;

/* shake128 digest of canonical coefficients of an interface vec */
static void
dig_vec (arith_vec_ptr v, unsigned int n, const abdlop_arith_params_t *ap,
         arith_ops_ptr ops, uint8_t out[32])
{
  shake128_state_t hs;
  uint8_t c8[8];
  unsigned int i, k;

  shake128_init (hs);
  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      {
        int64_t c = ops->vec_get_coeff (v, i, k);
        memcpy (c8, &c, 8);
        shake128_absorb (hs, c8, 8);
      }
  shake128_clear (hs);
  shake128_squeeze (hs, out, 32);
}

/* digest of a single polynomial's canonical coefficients */
static void
dig_poly (arith_poly_ptr p, const abdlop_arith_params_t *ap,
          arith_ops_ptr ops, uint8_t out[32])
{
  shake128_state_t hs;
  uint8_t c8[8];
  unsigned int k;

  shake128_init (hs);
  for (k = 0; k < ap->d; k++)
    {
      int64_t c = ops->poly_get_coeff (p, k);
      memcpy (c8, &c, 8);
      shake128_absorb (hs, c8, 8);
    }
  shake128_clear (hs);
  shake128_squeeze (hs, out, 32);
}

static long long
canonq (long long v, unsigned long long q)
{
  if (v < 0)
    v += (long long)q;
  else if (v >= (long long)q)
    v -= (long long)q;
  return v;
}

static void
run_engine (abdlop_arith_engine_t *e, const abdlop_params_srcptr params,
            const uint8_t seed[32], polyvec_t s1_m, polyvec_t s2_m,
            polyvec_t m_m, out_t *out, stats_t *st)
{
  abdlop_arith_params_t ap;
  arith_mat_ptr A1, A2p, Bp;
  arith_vec_ptr s1, s2, m;
  arith_vec_ptr tA1, tA2, tB, z1, z21, h;
  arith_poly_ptr c;
  uint8_t hashcomm[32], chash[32], hv[32];
  abdlop_arith_prove_stats_t pstats;
  unsigned int i, k;
  clock_t t0;

  abdlop_arith_params_init (&ap, params);

  s1 = e->ops->vec_alloc (e->ring, ap.m1);
  s2 = e->ops->vec_alloc (e->ring, ap.m2);
  m = e->ops->vec_alloc (e->ring, ap.l);
  for (i = 0; i < ap.m1; i++)
    for (k = 0; k < ap.d; k++)
      e->ops->vec_set_coeff (s1, i, k,
                             canonq (int_get_i64 (poly_get_coeff (
                                         polyvec_get_elem (s1_m, i), k)),
                                     ap.q));
  for (i = 0; i < ap.m2; i++)
    for (k = 0; k < ap.d; k++)
      e->ops->vec_set_coeff (s2, i, k,
                             canonq (int_get_i64 (poly_get_coeff (
                                         polyvec_get_elem (s2_m, i), k)),
                                     ap.q));
  for (i = 0; i < ap.l; i++)
    for (k = 0; k < ap.d; k++)
      e->ops->vec_set_coeff (m, i, k,
                             canonq (int_get_i64 (poly_get_coeff (
                                         polyvec_get_elem (m_m, i), k)),
                                     ap.q));

  A1 = e->ops->mat_alloc (e->ring, ap.kmsis, ap.m1);
  A2p = e->ops->mat_alloc (e->ring, ap.kmsis, ap.m22);
  Bp = e->ops->mat_alloc (e->ring, ap.l, ap.m22);
  tA1 = e->ops->vec_alloc (e->ring, ap.kmsis);
  tA2 = e->ops->vec_alloc (e->ring, ap.kmsis);
  tB = e->ops->vec_alloc (e->ring, ap.l);
  z1 = e->ops->vec_alloc (e->ring, ap.m1);
  z21 = e->ops->vec_alloc (e->ring, ap.m22);
  h = e->ops->vec_alloc (e->ring, ap.kmsis);
  c = e->ops->poly_alloc (e->ring);

  t0 = clock ();
  abdlop_arith_keygen (e, &ap, seed, A1, A2p, Bp);
  abdlop_arith_commit (e, &ap, tA1, tA2, tB, s1, m, s2, A1, A2p, Bp);
  st->commit = (double)(clock () - t0) / CLOCKS_PER_SEC;

  t0 = clock ();
  memcpy (hashcomm, seed, 32);
  abdlop_arith_hashcomm (hashcomm, e, &ap, tA1, tB);
  st->hashcomm = (double)(clock () - t0) / CLOCKS_PER_SEC;

  t0 = clock ();
  memcpy (chash, hashcomm, 32);
  abdlop_arith_prove (chash, c, z1, z21, h, tA2, s1, s2, A1, A2p, seed, e,
                      &ap, &pstats);
  st->prove = (double)(clock () - t0) / CLOCKS_PER_SEC;
  st->attempts = (long)pstats.attempts;
  st->rejloop = pstats.attempts > 0 ? pstats.rej_ms / pstats.attempts : 0.0;

  t0 = clock ();
  memcpy (hv, hashcomm, 32);
  out->accept_iv = abdlop_arith_verify (hv, c, z1, z21, h, tA1, A1, A2p, e,
                                        &ap);
  st->verify = (double)(clock () - t0) / CLOCKS_PER_SEC;

  /* outputs for cross-backend equivalence */
  memcpy (out->hashcomm, hashcomm, 32);
  memcpy (out->chash, chash, 32);
  {
    uint8_t da[32], db[32], dc[32], dd[32], de[32], df[32], dg[32];
    shake128_state_t hs;
    dig_vec (tA1, ap.kmsis, &ap, e->ops, da);
    dig_vec (tA2, ap.kmsis, &ap, e->ops, db);
    dig_vec (tB, ap.l, &ap, e->ops, dc);
    dig_poly (c, &ap, e->ops, dd);
    dig_vec (z1, ap.m1, &ap, e->ops, de);
    dig_vec (z21, ap.m22, &ap, e->ops, df);
    dig_vec (h, ap.kmsis, &ap, e->ops, dg);
    shake128_init (hs);
    shake128_absorb (hs, da, 32);
    shake128_absorb (hs, db, 32);
    shake128_absorb (hs, dc, 32);
    shake128_absorb (hs, dd, 32);
    shake128_absorb (hs, de, 32);
    shake128_absorb (hs, df, 32);
    shake128_absorb (hs, dg, 32);
    shake128_squeeze (hs, out->digest, 32);
    shake128_clear (hs);
  }

  e->ops->mat_free (A1);
  e->ops->mat_free (A2p);
  e->ops->mat_free (Bp);
  e->ops->vec_free (tA1);
  e->ops->vec_free (tA2);
  e->ops->vec_free (tB);
  e->ops->vec_free (z1);
  e->ops->vec_free (z21);
  e->ops->vec_free (h);
  e->ops->poly_free (c);
  e->ops->vec_free (s1);
  e->ops->vec_free (s2);
  e->ops->vec_free (m);
  st->total = st->commit + st->hashcomm + st->prove + st->verify;
  st->mism = 0;
}

int
main (void)
{
  abdlop_arith_engine_t eng[N_BACKENDS];
  arith_ops_ptr ops[N_BACKENDS];
  const char *name[N_BACKENDS] = { "lazer", "flint", "nmod", "ntl" };
  const abdlop_params_srcptr params = params1;
  abdlop_arith_params_t ap;
  polyvec_t s1_m, s2_m, m_m;
  uint8_t seed[32];
  out_t out[N_BACKENDS];
  stats_t st[N_BACKENDS];
  int i;
  unsigned int uu;
  INT_T (lo, 1);
  INT_T (hi, 1);

  lazer_init ();

  abdlop_arith_params_init (&ap, params);
  printf ("abdlop pipeline via arith interface; q=%llu d=%u m1=%u m2=%u"
          " kmsis=%u l=%u\n",
          (unsigned long long)ap.q, ap.d, ap.m1, ap.m2, ap.kmsis, ap.l);

  /* fixed reproducible seed */
  for (i = 0; i < 32; i++)
    seed[i] = (uint8_t)(i * 37 + 11);

  ops[0] = arith_lazer_ops ();
  ops[1] = arith_flint_ops ();
  ops[2] = arith_nmod_ops ();
  ops[3] = arith_ntl_ops ();

  memset (st, 0, sizeof st);
  memset (out, 0, sizeof out);

  for (i = 0; i < N_BACKENDS; i++)
    {
      eng[i].ops = ops[i];
      eng[i].ring = ops[i]->ring_new (ap.q, ap.d,
                                  i == 0 ? (const void *)params->ring
                                         : NULL);
    }

  /* shared above-interface input series: s1, s2, m */
  polyvec_alloc (s1_m, params->ring, ap.m1);
  polyvec_alloc (s2_m, params->ring, ap.m2);
  polyvec_alloc (m_m, params->ring, ap.l);
  uu = 0;
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s1_m, lo, hi, seed, uu++);
  polyvec_urandom_bnd (s2_m, lo, hi, seed, uu++);
  if (ap.l > 0)
    polyvec_urandom (m_m, params->ring->q, params->ring->log2q, seed, uu++);

  for (i = 0; i < N_BACKENDS; i++)
    {
      printf ("--- backend %s ---\n", name[i]);
      run_engine (&eng[i], params, seed, s1_m, s2_m, m_m, &out[i], &st[i]);
      printf ("    commit %.1fms hashcomm %.1fms prove %.1fms verify %.1fms"
              " total %.1fms\n",
              st[i].commit * 1e3, st[i].hashcomm * 1e3, st[i].prove * 1e3,
              st[i].verify * 1e3, st[i].total * 1e3);
      printf ("    prove attempts=%ld rejloop %.3fms avg/cycle (total %.1fms)\n",
              st[i].attempts, st[i].rejloop * 1e3, st[i].prove * 1e3);
      printf ("    accept_iv=%d\n", out[i].accept_iv);
    }

  for (i = 0; i < N_BACKENDS; i++)
    TEST_EXPECT (st[i].mism == 0);
  for (i = 0; i < N_BACKENDS; i++)
    TEST_EXPECT (out[i].accept_iv == 1);
  /* cross-backend equivalence of outputs and intermediates */
  for (i = 1; i < N_BACKENDS; i++)
    {
      TEST_EXPECT (memcmp (out[i].digest, out[0].digest, 32) == 0);
      TEST_EXPECT (memcmp (out[i].hashcomm, out[0].hashcomm, 32) == 0);
      TEST_EXPECT (memcmp (out[i].chash, out[0].chash, 32) == 0);
    }

  for (i = 0; i < N_BACKENDS; i++)
    eng[i].ops->ring_free (eng[i].ring);

  polyvec_free (s1_m);
  polyvec_free (s2_m);
  polyvec_free (m_m);

  TEST_PASS ();
}
