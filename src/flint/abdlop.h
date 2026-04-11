#pragma once

#include "../../lazer.h"
#include <flint/fq_nmod_types.h>

typedef struct {
  fmpz_t q;
  fmpz_t qminus1;
  fmpz_t m;
  fmpz_t mby2;
  fmpz_t gamma;
  fmpz_t gammaby2;
  fmpz_t pow2D;
  fmpz_t pow2Dby2;
  unsigned int D;
  int m_odd;
  unsigned int log2m;
} dcompress_params_flint_struct;
typedef dcompress_params_flint_struct dcompress_params_flint_t[1];

typedef struct {
  fq_nmod_ctx_struct *ring;
  dcompress_params_flint_t dcompress;
  /* dimensions  */
  unsigned int m1;   /* length of "short" message s1 */
  unsigned int m2;   /* length of randomness s2 */
  unsigned int l;    /* length of "large" message m */
  unsigned int lext; /* length of extension of m */
  unsigned int kmsis;
  /* norms */
  fmpz_t Bsqr; /* floor (B^2) */
  int64_t nu;      /* s2 uniform in [-nu,nu]*/
  int64_t omega;   /* challenges uniform in [-omega,omega], o(c)=c */
  unsigned int log2omega;
  uint64_t eta; /* sqrt(l1(o(c)*c)) <= eta XXX sqrt? */
  /* rejection sampling */
  int rej1;                /* do rejection sampling on s1 */
  unsigned int log2stdev1; /* stdev1 = 1.55 * 2^log2stdev1 */
  fmpz_t scM1;         /* scaled M1: round(M1 * 2^128) */
  fmpz_t stdev1sqr;
  int rej2;                /* do rejection sampling on s2 */
  unsigned int log2stdev2; /* stdev2 = 1.55 * 2^log2stdev2 */
  fmpz_t scM2;         /* scaled M2: round(M2 * 2^128) */
  fmpz_t stdev2sqr;
} abdlop_params_flint_struct;
typedef abdlop_params_flint_struct abdlop_params_flint_t[1];

void polymat_to_fq_mat(polymat_t src, fq_nmod_mat_t dst, fq_nmod_ctx_t ctx);
void polyvec_to_fq_mat(polyvec_t src, fq_nmod_mat_t dst, fq_nmod_ctx_t ctx);
void abdlop_params_to_flint(
    abdlop_params_flint_t dst,
    const abdlop_params_t src,
    fq_nmod_ctx_t ring
);

void fq_nmod_mat_addmul(fq_nmod_mat_t r, fq_nmod_mat_t a, fq_nmod_mat_t b, fq_nmod_ctx_t ctx);

void abdlop_commit_flint(
    fq_nmod_mat_t ta1,
    fq_nmod_mat_t ta2,
    fq_nmod_mat_t tb,
    fq_nmod_mat_t s1,
    fq_nmod_mat_t m,
    fq_nmod_mat_t s2,
    fq_nmod_mat_t a1,
    fq_nmod_mat_t a2prime,
    fq_nmod_mat_t bprime,
    fq_nmod_ctx_t ctx,
    abdlop_params_flint_t params
);
void abdlop_commit_flint2 (
    polyvec_t tA1,
    polyvec_t tA2,
    polyvec_t tB,
    polyvec_t s1,
    polyvec_t m,
    polyvec_t s2,
    polymat_t A1,
    polymat_t A2prime,
    polymat_t Bprime,
    const abdlop_params_t params
);
