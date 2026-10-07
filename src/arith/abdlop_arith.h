#ifndef ABDLOP_ARITH_H
#define ABDLOP_ARITH_H
/**
 * \file   abdlop_arith.h
 * \brief  abdlop keygen/commit/hashcomm/prove/verify over the pluggable
 *         arith interface.
 *
 * This module mirrors src/abdlop.c (the lazer pipeline), but every piece
 * of ring arithmetic runs through the ops of arith.h, so exactly the same
 * straight-line protocol produces identical values in every backend.
 *
 * Input generation (sampling, RNG, hashing) stays ABOVE this layer,
 * exactly as documented in arith.h: the public matrices are expanded from
 * the seed with the lazer uniform samplers and imported through the
 * coefficient accessors, and the prove loop drives the gaussian/autostable
 * samplers and the shake/coder machinery directly (they are deterministic,
 * backend independent, and only ever touch coefficient accessors).
 *
 * Everything compression/dcompression, rejection and norm related is done
 * coefficient by coefficient, in plain integer arithmetic, above the
 * interface; only multiplications and accumulating linear combinations go
 * through the backend ops.
 */

#include "lazer.h"
#include "arith.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* Per-backend engine: ops + concrete ring handle. */
  typedef struct
  {
    arith_ops_ptr ops;
    arith_ring_ptr ring;
  } abdlop_arith_engine_t;

  /* Scalar parameters extracted above the interface, so the pipeline is
   * backend independent. */
  typedef struct
  {
    unsigned long long q;
    unsigned int d, log2q;
    unsigned int D, log2m, log2qmd;
    unsigned int m1, m2, m22, kmsis, l;
    long long qminus1, m, mby2, gamma, gammaby2, pow2D, pow2Dby2, Bsqr;
    int m_odd;
    unsigned int omega, log2omega, log2stdev1, log2stdev2;
    int rej1, rej2;
    abdlop_params_srcptr p; /* lazer params: ring + scalars for sampling */
  } abdlop_arith_params_t;

  void abdlop_arith_params_init (abdlop_arith_params_t *ap,
                                 const abdlop_params_srcptr params);

  /* Expand uniformly random A1 (kmsis x m1), A2prime (kmsis x m22),
   * Bprime (l x m22) from seed||0, seed||1, seed||2 respectively.
   * Caller allocates A1, A2prime, Bprime with the arith ops of e. */
  void abdlop_arith_keygen (abdlop_arith_engine_t *e,
                            const abdlop_arith_params_t *ap,
                            const uint8_t seed[32], arith_mat_ptr A1,
                            arith_mat_ptr A2prime, arith_mat_ptr Bprime);

  /* tA1,tA2 in Rq^kmsis, tB in Rq^l.
   * tA2 = s22 + A1*s1 + A2prime*s21 ;
   * tA1 = Power2Round(tA2,D) ; tA2 <- tA2 - 2^D*tA1 ;
   * tB  = m + Bprime*s21.
   * Caller allocates tA1,tA2,tB and imports s1 (m1), m (l), s2 (m2),
   * A1, A2prime, Bprime through the accessors. */
  void abdlop_arith_commit (abdlop_arith_engine_t *e,
                            const abdlop_arith_params_t *ap,
                            arith_vec_ptr tA1, arith_vec_ptr tA2,
                            arith_vec_ptr tB, const arith_vec_ptr s1,
                            const arith_vec_ptr m, const arith_vec_ptr s2,
                            const arith_mat_ptr A1,
                            const arith_mat_ptr A2prime,
                            const arith_mat_ptr Bprime);

  /* hash <- H(hash || encode(tA1,tB)), matching abdlop_hashcomm. */
  void abdlop_arith_hashcomm (uint8_t hash[32], abdlop_arith_engine_t *e,
                              const abdlop_arith_params_t *ap,
                              const arith_vec_ptr tA1,
                              const arith_vec_ptr tB);

  /* Compress an opening proof (c,z1,z21,h) from hash of transcript, short
   * message s1 and randomness s2, tA2 and the public key parts.
   * Rejection sampling is deterministic (see params) and driven with the
   * lazer samplers above the interface; updates hash to the challenge seed.
   * Caller allocates c, z1 (m1), z21 (m22), h (kmsis) and imports tA2,
   * s1, s2, A1, A2prime.  stats, if non-NULL, receives the rejection loop
   * breakdown: attempts counts the rejection cycles (loop iterations), rej_ms
   * is the time inside the loop (sampling + arithmetic + rejection) and
   * post_ms the time after the accepted attempt (w1post/hint/finalize). */
  typedef struct
  {
    unsigned int attempts;
    double rej_ms;
    double post_ms;
  } abdlop_arith_prove_stats_t;

  void abdlop_arith_prove (uint8_t hash[32], arith_poly_ptr c,
                           arith_vec_ptr z1, arith_vec_ptr z21,
                           arith_vec_ptr h, const arith_vec_ptr tA2,
                           const arith_vec_ptr s1, const arith_vec_ptr s2,
                           const arith_mat_ptr A1,
                           const arith_mat_ptr A2prime,
                           const uint8_t seed[32],
                           abdlop_arith_engine_t *e,
                           const abdlop_arith_params_t *ap,
                           abdlop_arith_prove_stats_t *stats);

  /* Verify an opening proof (c,z1,z21,h) from hash of transcript,
   * commitment part tA1 and public key parts; updates hash.
   * Returns 1 on successful verification. */
  int abdlop_arith_verify (uint8_t hash[32], const arith_poly_ptr c,
                           const arith_vec_ptr z1, const arith_vec_ptr z21,
                           const arith_vec_ptr h, const arith_vec_ptr tA1,
                           const arith_mat_ptr A1,
                           const arith_mat_ptr A2prime,
                           abdlop_arith_engine_t *e,
                           const abdlop_arith_params_t *ap);

#ifdef __cplusplus
}
#endif

#endif /* ABDLOP_ARITH_H */