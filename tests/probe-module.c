/*
 * \file   probe-module.c
 * \brief  Compare abdlop_arith_* (via the lazer interface backend) against
 *         the lazer reference pipeline on the SAME seed/inputs.
 *
 * One-off validation probe: every intermediate and output of the module
 * is cross-checked against tests/abdlop-test.c semantics.  Not part of
 * the test suite.
 */

#include "abdlop-params1.h"
#include "lazer.h"

#include "../src/arith/arith.h"
#include "../src/arith/abdlop_arith.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "test.h"

static long long
canonq (long long v, unsigned long long q)
{
  if (v < 0)
    v += (long long)q;
  else if (v >= (long long)q)
    v -= (long long)q;
  return v;
}

int
main (void)
{
  const abdlop_params_srcptr params = params1;
  polyring_srcptr Rq = params->ring;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int m2 = params->m2;
  const unsigned int m22 = m2 - kmsis;
  const unsigned int l = params->l + params->lext;
  const unsigned int d = Rq->d;
  uint8_t seed[32];
  int i, k;
  unsigned int uu;
  long mism = 0;

  for (i = 0; i < 32; i++)
    seed[i] = (uint8_t)(i * 37 + 11);

  lazer_init ();

  /* ---------------- reference (lazer) objects ---------------- */
  poly_t c_ref, c2_ref;
  polyvec_t s1_r, s2_r, m_r, tA1_r, tA2_r, tB_r;
  polyvec_t z1_r, z21_r, h_r;
  polymat_t A1_r, A2p_r, Bp_r;
  uint8_t hashp[32], hashv[32], hashc[32];
  INT_T (lo, 1);
  INT_T (hi, 1);

  poly_alloc (c_ref, Rq);
  poly_alloc (c2_ref, Rq);
  polyvec_alloc (s1_r, Rq, m1);
  polyvec_alloc (s2_r, Rq, m2);
  polyvec_alloc (m_r, Rq, l);
  polyvec_alloc (tA1_r, Rq, kmsis);
  polyvec_alloc (tA2_r, Rq, kmsis);
  polyvec_alloc (tB_r, Rq, l);
  polyvec_alloc (z1_r, Rq, m1);
  polyvec_alloc (z21_r, Rq, m22);
  polyvec_alloc (h_r, Rq, kmsis);
  polymat_alloc (A1_r, Rq, kmsis, m1);
  polymat_alloc (A2p_r, Rq, kmsis, m22);
  polymat_alloc (Bp_r, Rq, l, m22);

  uu = 0;
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s1_r, lo, hi, seed, uu++);
  polyvec_urandom_bnd (s2_r, lo, hi, seed, uu++);
  if (l > 0)
    polyvec_urandom (m_r, Rq->q, Rq->log2q, seed, uu++);

  abdlop_keygen (A1_r, A2p_r, Bp_r, seed, params);
  abdlop_commit (tA1_r, tA2_r, tB_r, s1_r, m_r, s2_r, A1_r, A2p_r, Bp_r,
                 params);

  memcpy (hashc, seed, 32);
  abdlop_hashcomm (hashc, tA1_r, tB_r, params);
  memcpy (hashp, hashc, 32);
  abdlop_prove (hashp, c_ref, z1_r, z21_r, h_r, tA2_r, s1_r, s2_r, A1_r,
                A2p_r, seed, params);
  memcpy (hashv, seed, 32);
  abdlop_hashcomm (hashv, tA1_r, tB_r, params);
  {
    int b = abdlop_verify (hashv, c_ref, z1_r, z21_r, h_r, tA1_r, A1_r,
                           A2p_r, params);
    printf ("ref verify accept=%d hash_eq=%d\n", b,
            (int)(memcmp (hashp, hashv, 32) == 0));
  }

  /* ---------------- module (lazer interface backend) ---------------- */
  abdlop_arith_engine_t eng;
  abdlop_arith_params_t ap;
  arith_vec_ptr s1_i, s2_i, m_i, tA1_i, tA2_i, tB_i, z1_i, z21_i, h_i;
  arith_mat_ptr A1_i, A2p_i, Bp_i;
  arith_poly_ptr c_i;
  uint8_t hash_p[32], hash_v[32];

  abdlop_arith_params_init (&ap, params);
  eng.ops = arith_lazer_ops ();
  eng.ring = eng.ops->ring_new (ap.q, ap.d, (const void *)Rq);

  s1_i = eng.ops->vec_alloc (eng.ring, m1);
  s2_i = eng.ops->vec_alloc (eng.ring, m2);
  m_i = eng.ops->vec_alloc (eng.ring, l);
  tA1_i = eng.ops->vec_alloc (eng.ring, kmsis);
  tA2_i = eng.ops->vec_alloc (eng.ring, kmsis);
  tB_i = eng.ops->vec_alloc (eng.ring, l);
  z1_i = eng.ops->vec_alloc (eng.ring, m1);
  z21_i = eng.ops->vec_alloc (eng.ring, m22);
  h_i = eng.ops->vec_alloc (eng.ring, kmsis);
  A1_i = eng.ops->mat_alloc (eng.ring, kmsis, m1);
  A2p_i = eng.ops->mat_alloc (eng.ring, kmsis, m22);
  Bp_i = eng.ops->mat_alloc (eng.ring, l, m22);
  c_i = eng.ops->poly_alloc (eng.ring);

  /* copy reference s1/s2/m into the interface */
  for (i = 0; i < (int)m1; i++)
    for (k = 0; k < (int)d; k++)
      eng.ops->vec_set_coeff (s1_i, i, k,
                              canonq (int_get_i64 (poly_get_coeff (
                                          polyvec_get_elem (s1_r, i), k)),
                                      ap.q));
  for (i = 0; i < (int)m2; i++)
    for (k = 0; k < (int)d; k++)
      eng.ops->vec_set_coeff (s2_i, i, k,
                              canonq (int_get_i64 (poly_get_coeff (
                                          polyvec_get_elem (s2_r, i), k)),
                                      ap.q));
  for (i = 0; i < (int)l; i++)
    for (k = 0; k < (int)d; k++)
      eng.ops->vec_set_coeff (m_i, i, k,
                              canonq (int_get_i64 (poly_get_coeff (
                                          polyvec_get_elem (m_r, i), k)),
                                      ap.q));

  /* module keygen: must reproduce the reference matrices */
  abdlop_arith_keygen (&eng, &ap, seed, A1_i, A2p_i, Bp_i);
  for (i = 0; i < (int)kmsis; i++)
    for (k = 0; k < (int)m1; k++)
      for (uu = 0; uu < d; uu++)
        {
          long long a = eng.ops->mat_get_coeff (A1_i, i, k, uu);
          long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polymat_get_elem (A1_r, i, k), uu)),
                                ap.q);
          if (a != b)
            {
              if (mism < 5)
                printf ("    [kg-A1] %d %d %u iv %lld ref %lld\n", i, k, uu,
                        (long long)a, (long long)b);
              mism++;
            }
        }
  for (i = 0; i < (int)kmsis; i++)
    for (k = 0; k < (int)m22; k++)
      for (uu = 0; uu < d; uu++)
        {
          long long a = eng.ops->mat_get_coeff (A2p_i, i, k, uu);
          long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polymat_get_elem (A2p_r, i, k), uu)),
                                ap.q);
          if (a != b)
            {
              if (mism < 5)
                printf ("    [kg-A2p] %d %d %u iv %lld ref %lld\n", i, k, uu,
                        (long long)a, (long long)b);
              mism++;
            }
        }
  for (i = 0; i < (int)l; i++)
    for (k = 0; k < (int)m22; k++)
      for (uu = 0; uu < d; uu++)
        {
          long long a = eng.ops->mat_get_coeff (Bp_i, i, k, uu);
          long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polymat_get_elem (Bp_r, i, k), uu)),
                                ap.q);
          if (a != b)
            {
              if (mism < 5)
                printf ("    [kg-Bp] %d %d %u iv %lld ref %lld\n", i, k, uu,
                        (long long)a, (long long)b);
              mism++;
            }
        }

  /* module commit */
  abdlop_arith_commit (&eng, &ap, tA1_i, tA2_i, tB_i, s1_i, m_i, s2_i,
                       A1_i, A2p_i, Bp_i);
  for (i = 0; i < (int)kmsis; i++)
    for (k = 0; k < (int)d; k++)
      {
        long long a = eng.ops->vec_get_coeff (tA1_i, i, k);
        long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polyvec_get_elem (tA1_r, i), k)),
                              ap.q);
        if (a != b)
          {
            if (mism < 5)
              printf ("    [cm-tA1] %d %d iv %lld ref %lld\n", i, k,
                      (long long)a, (long long)b);
            mism++;
          }
      }
  for (i = 0; i < (int)kmsis; i++)
    for (k = 0; k < (int)d; k++)
      {
        long long a = eng.ops->vec_get_coeff (tA2_i, i, k);
        long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polyvec_get_elem (tA2_r, i), k)),
                              ap.q);
        if (a != b)
          {
            if (mism < 5)
              printf ("    [cm-tA2] %d %d iv %lld ref %lld\n", i, k,
                      (long long)a, (long long)b);
            mism++;
          }
      }
  for (i = 0; i < (int)l; i++)
    for (k = 0; k < (int)d; k++)
      {
        long long a = eng.ops->vec_get_coeff (tB_i, i, k);
        long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polyvec_get_elem (tB_r, i), k)),
                              ap.q);
        if (a != b)
          {
            if (mism < 5)
              printf ("    [cm-tB] %d %d iv %lld ref %lld\n", i, k,
                      (long long)a, (long long)b);
            mism++;
          }
      }

  /* module hashcomm */
  memcpy (hash_p, seed, 32);
  abdlop_arith_hashcomm (hash_p, &eng, &ap, tA1_i, tB_i);
  if (memcmp (hash_p, hashc, 32) != 0)
    {
      printf ("    [hashcomm] module != ref\n");
      mism++;
    }

  /* module prove */
  abdlop_arith_prove (hash_p, c_i, z1_i, z21_i, h_i, tA2_i, s1_i, s2_i,
                      A1_i, A2p_i, seed, &eng, &ap, NULL);
  if (memcmp (hash_p, hashp, 32) != 0)
    {
      printf ("    [prove-hash] module != ref\n");
      mism++;
    }
  for (k = 0; k < (int)d; k++)
    {
      long long a = eng.ops->poly_get_coeff (c_i, k);
      long long b = canonq (int_get_i64 (poly_get_coeff (c_ref, k)), ap.q);
      if (a != b)
        {
          if (mism < 5)
            printf ("    [prove-c] k %d iv %lld ref %lld\n", k, (long long)a,
                    (long long)b);
          mism++;
        }
    }
  for (i = 0; i < (int)m1; i++)
    for (k = 0; k < (int)d; k++)
      {
        long long a = eng.ops->vec_get_coeff (z1_i, i, k);
        long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polyvec_get_elem (z1_r, i), k)),
                              ap.q);
        if (a != b)
          {
            if (mism < 5)
              printf ("    [prove-z1] %d %d iv %lld ref %lld\n", i, k,
                      (long long)a, (long long)b);
            mism++;
          }
      }
  for (i = 0; i < (int)m22; i++)
    for (k = 0; k < (int)d; k++)
      {
        long long a = eng.ops->vec_get_coeff (z21_i, i, k);
        long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polyvec_get_elem (z21_r, i), k)),
                              ap.q);
        if (a != b)
          {
            if (mism < 5)
              printf ("    [prove-z21] %d %d iv %lld ref %lld\n", i, k,
                      (long long)a, (long long)b);
            mism++;
          }
      }
  for (i = 0; i < (int)kmsis; i++)
    for (k = 0; k < (int)d; k++)
      {
        long long a = eng.ops->vec_get_coeff (h_i, i, k);
        long long b = canonq (int_get_i64 (poly_get_coeff (
                                    polyvec_get_elem (h_r, i), k)),
                              ap.q);
        if (a != b)
          {
            if (mism < 5)
              printf ("    [prove-h] %d %d iv %lld ref %lld\n", i, k,
                      (long long)a, (long long)b);
            mism++;
          }
      }

  /* module verify */
  memcpy (hash_v, seed, 32);
  abdlop_arith_hashcomm (hash_v, &eng, &ap, tA1_i, tB_i);
  {
    int b = abdlop_arith_verify (hash_v, c_i, z1_i, z21_i, h_i, tA1_i,
                                 A1_i, A2p_i, &eng, &ap);
    printf ("module verify accept=%d hash_eq=%d\n", b,
            (int)(memcmp (hash_v, hashp, 32) == 0));
  }

  printf ("module-vs-ref mism=%ld\n", mism);

  /* cleanup */
  eng.ops->vec_free (s1_i);
  eng.ops->vec_free (s2_i);
  eng.ops->vec_free (m_i);
  eng.ops->vec_free (tA1_i);
  eng.ops->vec_free (tA2_i);
  eng.ops->vec_free (tB_i);
  eng.ops->vec_free (z1_i);
  eng.ops->vec_free (z21_i);
  eng.ops->vec_free (h_i);
  eng.ops->mat_free (A1_i);
  eng.ops->mat_free (A2p_i);
  eng.ops->mat_free (Bp_i);
  eng.ops->poly_free (c_i);
  eng.ops->ring_free (eng.ring);

  poly_free (c_ref);
  poly_free (c2_ref);
  polyvec_free (s1_r);
  polyvec_free (s2_r);
  polyvec_free (m_r);
  polyvec_free (tA1_r);
  polyvec_free (tA2_r);
  polyvec_free (tB_r);
  polyvec_free (z1_r);
  polyvec_free (z21_r);
  polyvec_free (h_r);
  polymat_free (A1_r);
  polymat_free (A2p_r);
  polymat_free (Bp_r);

  if (mism == 0)
    printf ("PROBE OK\n");
  return mism == 0 ? 0 : 1;
}