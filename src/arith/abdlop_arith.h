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
 * Input generation (sampling, RNG, hashing) is handled internally with
 * self-contained, deterministic implementations (no lazer.h dependency);
 * only the arith interface is required.
 *
 * Everything compression/dcompression, rejection and norm related is done
 * coefficient by coefficient, in plain integer arithmetic, above the
 * interface; only multiplications and accumulating linear combinations go
 * through the backend ops.
 */

#include "arith.h"
#include <stdint.h>

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

  /* Scalar parameters, fully self-contained (no lazer types). */
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
    /* rejection sampling: scM = round(M * 2^128) as 3 limbs (le), sigma^2 */
    uint64_t scM1[3], scM2[3];
    uint64_t stdev1sqr[2], stdev2sqr[2]; /* 128-bit enough, nlimbs<=2 */
  } abdlop_arith_params_t;

  /* Initialize scalar params from explicit values. */
  void abdlop_arith_params_set (abdlop_arith_params_t *ap,
                                unsigned long long q, unsigned int d,
                                unsigned int log2q, unsigned int D,
                                unsigned int log2m, unsigned int m1,
                                unsigned int m2, unsigned int kmsis,
                                unsigned int l, unsigned int lext,
                                long long qminus1, long long mval,
                                long long mby2v, long long gammav,
                                long long gammaby2v, long long pow2Dv,
                                long long pow2Dby2v, long long Bsqr,
                                int m_odd, unsigned int omega,
                                unsigned int log2omega,
                                unsigned int log2stdev1, int rej1,
                                const uint64_t scM1[3],
                                const uint64_t stdev1sqr[2],
                                unsigned int log2stdev2, int rej2,
                                const uint64_t scM2[3],
                                const uint64_t stdev2sqr[2]);

  /* Expand uniformly random A1 (kmsis x m1), A2prime (kmsis x m22),
   * Bprime (l x m22) from seed||0, seed||1, seed||2 respectively.
   * Caller allocates A1, A2prime, Bprime with the arith ops of e. */
  void abdlop_arith_keygen (abdlop_arith_engine_t *e,
                            const abdlop_arith_params_t *ap,
                            const uint8_t seed[32], arith_mat_ptr A1,
                            arith_mat_ptr A2prime, arith_mat_ptr Bprime);

  /* tA1,tA2 in Rq^kmsis, tB in Rq^l. */
  void abdlop_arith_commit (abdlop_arith_engine_t *e,
                            const abdlop_arith_params_t *ap,
                            arith_vec_ptr tA1, arith_vec_ptr tA2,
                            arith_vec_ptr tB, const arith_vec_ptr s1,
                            const arith_vec_ptr m, const arith_vec_ptr s2,
                            const arith_mat_ptr A1,
                            const arith_mat_ptr A2prime,
                            const arith_mat_ptr Bprime);

  /* hash <- H(hash || encode(tA1,tB)). */
  void abdlop_arith_hashcomm (uint8_t hash[32], abdlop_arith_engine_t *e,
                              const abdlop_arith_params_t *ap,
                              const arith_vec_ptr tA1,
                              const arith_vec_ptr tB);

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
