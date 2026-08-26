#include "abdlop.h"
#include "src/dom.h"
#include "src/flint/abdlop_coder.h"
#include "src/intvec.h"
#include "src/urandom.h"
#include <flint/fmpz_extras.h>
#include <mpfr.h>
#include <flint/flint.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <time.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_types.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/arb.h>

void abdlop_commit_flint(
    fq_default_mat_t tA1,
    fq_default_mat_t tA2,
    fq_default_mat_t tB,
    fq_default_mat_t s1,
    fq_default_mat_t m,
    fq_default_mat_t s2,
    fq_default_mat_t A1,
    fq_default_mat_t A2prime,
    fq_default_mat_t Bprime,
    const abdlop_params_flint_t params
) {
    // Extract parameters
    fq_default_ctx_struct *ctx = params->ring;
    const unsigned int kmsis = params->kmsis;
    const unsigned int m1 = params->m1;
    const unsigned int m2 = params->m2;
    const unsigned int l = params->l;

    // 1. Setup Windows for s2
    fq_default_mat_t s21, s22;
    // fq_default_mat_window_init(window, source, r1, c1, r2, c2, ctx)
    fq_default_mat_window_init(s21, s2, 0, 0, m2 - kmsis, 1, ctx);
    fq_default_mat_window_init(s22, s2, m2 - kmsis, 0, m2, 1, ctx);

    if (m1 > 0) {
        // 2. Compute tA2 = s22 + A1*s1 + A2prime*s21
        // We use tA2 as the accumulator directly
        fq_default_mat_set(tA2, s22, ctx);

        fq_default_mat_t tmp_kmsis;
        fq_default_mat_init(tmp_kmsis, kmsis, 1, ctx);

        // tA2 = tA2 + A1 * s1
        fq_default_mat_mul(tmp_kmsis, A1, s1, ctx);
        fq_default_mat_add(tA2, tA2, tmp_kmsis, ctx);

        // tA2 = tA2 + A2prime * s21
        fq_default_mat_mul(tmp_kmsis, A2prime, s21, ctx);
        fq_default_mat_add(tA2, tA2, tmp_kmsis, ctx);

        // printf("\n"); fq_default_mat_print_pretty(tA2, ctx); printf("\n");
        fq_default_mat_power2round(tA1, tA2, params);
        // printf("\n"); fq_default_mat_print_pretty(tA2, ctx); printf("\n");
        // printf("\n"); fq_default_mat_print_pretty(tA1, ctx); printf("\n");

        fq_default_mat_clear(tmp_kmsis, ctx);
    }
    if (l > 0) {
        // 4. Handle Message Commitment tB
        fq_default_mat_t m_sub, tB_sub, Bprime_sub, tmp_l;

        // Initialize windows for sub-sections of m, tB, and Bprime
        fq_default_mat_window_init(m_sub, m, 0, 0, l, 1, ctx);
        fq_default_mat_window_init(tB_sub, tB, 0, 0, l, 1, ctx);
        fq_default_mat_window_init(Bprime_sub, Bprime, 0, 0, l, m2 - kmsis, ctx);

        fq_default_mat_init(tmp_l, l, 1, ctx);

        // tB_sub = Bprime_sub * s21 + m_sub
        fq_default_mat_mul(tmp_l, Bprime_sub, s21, ctx);
        fq_default_mat_add(tB_sub, tmp_l, m_sub, ctx);

        // Cleanup windows and temp for this block
        fq_default_mat_clear(tmp_l, ctx);
        fq_default_mat_window_clear(m_sub, ctx);
        fq_default_mat_window_clear(tB_sub, ctx);
        fq_default_mat_window_clear(Bprime_sub, ctx);
    }

    // Final Cleanup
    fq_default_mat_window_clear(s21, ctx);
    fq_default_mat_window_clear(s22, ctx);
}

void abdlop_enccomm_flint(
    uint8_t *buf,
    size_t *buflen,
    fq_default_mat_t tA1,
    fq_default_mat_t tB,
    abdlop_params_flint_t params
) {
  slong log2q = fmpz_clog_ui(params->dcompress->qminus1, 2);
  const unsigned int d = fq_default_ctx_degree(params->ring);
  const unsigned int D = params->dcompress->D;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int l = params->l;

  coder_state_t cstate;
  fq_default_mat_t tb_;
  const unsigned int len
      = CEIL (kmsis * d * (log2q - D) + l * d * log2q, 8) + 1;

  if (buflen != NULL)
    *buflen = len;

  if (buf == NULL)
    return;

  coder_enc_begin (cstate, buf);
  if (m1 > 0) {
      fmpz_t mod;
      fmpz_init(mod);
      fq_default_ctx_prime(mod, params->ring);

      fmpz_set_ui(mod, 1);
      fmpz_mul_2exp(mod, mod, log2q - D);

      coder_enc_urandom3_flint(
          cstate,
          tA1,
          params->ring,
          params->mod_ctx,
          params->dcompress->m,
          log2q - D
      );
  }

  if (l > 0) {
      fq_default_mat_window_init(tb_, tB, 0, 0, l, 1, params->ring);
      coder_enc_urandom3_flint(
          cstate,
          tb_,
          params->ring,
          params->mod_ctx,
          params->dcompress->m,
          log2q
      );
  }

  coder_enc_end (cstate);
}

void abdlop_hashcomm_flint(
    uint8_t hash[32],
    fq_default_mat_t ta1,
    fq_default_mat_t tb,
    abdlop_params_flint_t params
) {
    slong log2q = fmpz_clog_ui(params->dcompress->qminus1, 2);
    const unsigned int d = fq_default_ctx_degree(params->ring);
    const unsigned int D = params->dcompress->D;
    const unsigned int kmsis = params->kmsis;
    const unsigned int l = params->l;
    shake128_state_t hstate;
    const size_t outlen = CEIL (kmsis * d * (log2q - D) + l * d * log2q, 8) + 1;
    uint8_t out[outlen];

    abdlop_enccomm_flint(
        out,
        NULL,
        ta1,
        tb,
        params
    );

    shake128_init (hstate);
    shake128_absorb (hstate, hash, 32);
    shake128_absorb (hstate, out, outlen);
    shake128_squeeze (hstate, hash, 32);
    shake128_clear (hstate);
}

void abdlop_prove_flint(
    uint8_t hash[32],
    fmpz_mod_poly_t c,
    fq_default_mat_t z1,
    fq_default_mat_t z21,
    fq_default_mat_t h,
    fq_default_mat_t tA2,
    fq_default_mat_t s1,
    fq_default_mat_t s2,
    fq_default_mat_t A1,
    fq_default_mat_t A2prime,
    const uint8_t seed[32],
    abdlop_params_flint_t params
) {
    slong kmsis = params->kmsis;
    slong m1 = params->m1;
    slong m2 = params->m2;
    slong log2m = params->dcompress->log2m;

    fmpz_mod_poly_t modulus_poly;
    fmpz_mod_poly_init(modulus_poly, params->mod_ctx);
    fq_default_ctx_modulus(modulus_poly, params->ring);
    slong d = fmpz_mod_poly_length(modulus_poly, params->mod_ctx) - 1;
    fmpz_mod_poly_clear(modulus_poly, params->mod_ctx);

    fq_default_mat_t y1, y2, cs1, cs2, w, w1, w0;
    fq_default_mat_init(y1, m1, 1, params->ring);
    fq_default_mat_init(y2, m2, 1, params->ring);
    fq_default_mat_init(cs1, m1, 1, params->ring);
    fq_default_mat_init(cs2, m2, 1, params->ring);
    fq_default_mat_init(w, kmsis, 1, params->ring);
    fq_default_mat_init(w1, kmsis, 1, params->ring);
    fq_default_mat_init(w0, kmsis, 1, params->ring);

    fq_default_mat_t s21, s22;
    fq_default_mat_t y21, y22;

    // Correct Window initializations over parent matrices
    fq_default_mat_window_init(s21, s2, 0, 0, m2 - kmsis, 1, params->ring);
    fq_default_mat_window_init(s22, s2, m2 - kmsis, 0, m2, 1, params->ring);

    fq_default_mat_window_init(y21, y2, 0, 0, m2 - kmsis, 1, params->ring);
    fq_default_mat_window_init(y22, y2, m2 - kmsis, 0, m2, 1, params->ring);

    uint8_t cseed[32];
    uint8_t yseed[32];
    rng_state_t rngstate;
    rng_init(rngstate, seed, 0);
    rng_urandom(rngstate, yseed, 32);

    coder_state_t cstate;
    unsigned int outlen;
    uint8_t out[CEIL(log2m * d * kmsis, 8) + 1];

    shake128_state_t hstate;

    int rej;

    fq_default_t cd;
    fq_default_init2(cd, params->ring);

    fq_default_mat_t tmp;
    fq_default_mat_init(tmp, kmsis, 1, params->ring);

    fq_default_mat_t tmp_mul;
    fq_default_mat_init(tmp_mul, kmsis, 1, params->ring);

    uint32_t dom = 0;
    while (1) {
        fq_default_mat_grand(y1, params->ring, params->mod_ctx, params->dcompress->q, params->log2stdev1, yseed, dom);
        dom++;
        fq_default_mat_grand(y2, params->ring, params->mod_ctx, params->dcompress->q, params->log2stdev2, yseed, dom);
        dom++;

        // w = y22 + A1*y1 + A2prime*y21
        fq_default_mat_set(w, y22, params->ring);

        fq_default_mat_mul(tmp_mul, A1, y1, params->ring);
        fq_default_mat_add(w, w, tmp_mul, params->ring);

        fq_default_mat_mul(tmp_mul, A2prime, y21, params->ring);
        fq_default_mat_add(w, w, tmp_mul, params->ring);

        fq_default_mat_dcompress_decompose(w1, w0, w, params);

        coder_enc_begin(cstate, out);
        coder_enc_urandom3_flint(
            cstate,
            w1,
            params->ring,
            params->mod_ctx,
            params->dcompress->m,
            log2m
        );
        coder_enc_end(cstate);

        outlen = coder_get_offset(cstate);
        ASSERT_ERR(outlen % 8 == 0);
        ASSERT_ERR(outlen / 8 <= CEIL(log2m * d * kmsis, 8) + 1);
        outlen >>= 3; /* nbits to nbytes */

        shake128_init(hstate);
        shake128_absorb(hstate, hash, 32);
        shake128_absorb(hstate, out, outlen);
        shake128_squeeze(hstate, cseed, 32);

        fmpz_mod_poly_urand_autostable(c, d, params->mod_ctx, params->omega, params->log2omega, cseed, 0);

        fq_default_set_fmpz_mod_poly(cd, c, params->ring);
        fq_default_mat_scalar_mul(cs1, s1, cd, params->ring);
        fq_default_mat_scalar_mul(cs2, s2, cd, params->ring);
        fq_default_mat_add(y1, y1, cs1, params->ring);
        fq_default_mat_add(y2, y2, cs2, params->ring);

        if (params->rej1) {
            rej = rej_standard_flint(
                rngstate,
                y1,
                cs1,
                params->ring,
                params->mod_ctx,
                params->scM1,
                params->stdev1sqr
            );
            if (rej) continue;
        }
        if (params->rej2) {
            rej = rej_standard_flint(
                rngstate,
                y2,
                cs2,
                params->ring,
                params->mod_ctx,
                params->scM2,
                params->stdev2sqr
            );
            if (rej) continue;
        }

        // y22 = y22 - c*tA2 - w0
        fq_default_mat_scalar_mul(tmp, tA2, cd, params->ring);
        fq_default_mat_sub(y22, y22, tmp, params->ring);
        fq_default_mat_sub(y22, y22, w0, params->ring);

        fmpz_t norm;
        fmpz_init(norm);
        polyvec_l2sqr_flint(norm, y2, params->ring, params->mod_ctx);
        rej = fmpz_cmp(norm, params->Bsqr) > 0;
        fmpz_clear(norm);
        if (rej) continue;

        break;
    }

    fmpz_mod_poly_t gamma_poly;
    fmpz_mod_poly_init(gamma_poly, params->mod_ctx);
    fmpz_mod_poly_set_fmpz(gamma_poly, params->dcompress->gamma, params->mod_ctx);

    fq_default_t g;
    fq_default_init2(g, params->ring);
    fq_default_set_fmpz_mod_poly(g, gamma_poly, params->ring);
    fmpz_mod_poly_clear(gamma_poly, params->mod_ctx);

    // w1 = gamma * w1 - y22
    fq_default_mat_scalar_mul(w1, w1, g, params->ring);
    fq_default_clear(g, params->ring);

    fq_default_mat_sub(w1, w1, y22, params->ring);

    // Write results to output buffers
    fq_default_mat_set(z1, y1, params->ring);
    fq_default_mat_set(z21, y21, params->ring);

    fq_default_mat_dcompress_make_ghint(h, y22, w1, params);

    memcpy(hash, cseed, 32);

    /* Clean up resources */
    shake128_clear(hstate);
    rng_clear(rngstate);
    fq_default_mat_clear(y1, params->ring);
    fq_default_mat_clear(y2, params->ring);
    fq_default_mat_clear(cs1, params->ring);
    fq_default_mat_clear(cs2, params->ring);
    fq_default_mat_clear(w, params->ring);
    fq_default_mat_clear(w1, params->ring);
    fq_default_mat_clear(w0, params->ring);
    fq_default_mat_clear(tmp, params->ring);
    fq_default_mat_clear(tmp_mul, params->ring);
    fq_default_clear(cd, params->ring);

    fq_default_mat_window_clear(s21, params->ring);
    fq_default_mat_window_clear(s22, params->ring);
    fq_default_mat_window_clear(y21, params->ring);
    fq_default_mat_window_clear(y22, params->ring);
}

void abdlop_keygen_flint(
    fq_default_mat_t A1,
    fq_default_mat_t A2prime,
    fq_default_mat_t Bprime,
    const uint8_t seed[32],
    abdlop_params_flint_t params
) {
    slong log2q = fmpz_clog_ui(params->dcompress->qminus1, 2);
    if (params->m1 > 0) {
        fq_default_mat_urand(
            A1,
            params->ring,
            params->mod_ctx,
            params->dcompress->q,
            log2q,
            seed,
            0
        );
        fq_default_mat_urand(
            A2prime,
            params->ring,
            params->mod_ctx,
            params->dcompress->q,
            log2q,
            seed,
            1
        );
        // printf("\n"); fq_default_mat_print_pretty(A1, params->ring); printf("\n");
        // printf("\n"); fq_default_mat_print_pretty(A2prime, params->ring); printf("\n");
    }
    if (params->l + params->lext > 0) {
        fq_default_mat_urand(
            Bprime,
            params->ring,
            params->mod_ctx,
            params->dcompress->q,
            log2q,
            seed,
            2
        );
        // printf("\n"); fq_default_mat_print_pretty(Bprime, params->ring); printf("\n");
    }
}

void abdlop_keygen_flint2(
    polymat_t A1,
    polymat_t A2prime,
    polymat_t Bprime,
    const uint8_t seed[32],
    const abdlop_params_t params
) {
    fq_default_mat_t a1, a2prime, bprime;
    fmpz_mod_ctx_t mod_ctx;
    fq_default_ctx_t ctx;
    abdlop_params_flint_t fparams;
    fq_default_mat_init_polymat(A1, a1, ctx, mod_ctx);
    fq_default_mat_init_polymat(A2prime, a2prime, ctx, mod_ctx);
    fq_default_mat_init_polymat(Bprime, bprime, ctx, mod_ctx);
    abdlop_params_to_flint(fparams, params, ctx, mod_ctx);
    abdlop_keygen_flint(a1, a2prime, bprime, seed, fparams);
}

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
) {
    fq_default_mat_t ta1, ta2, tb, ms1, mm, ms2, a1, a2prime, bprime;
    fmpz_mod_ctx_t mod_ctx;
    fq_default_ctx_t ctx;
    abdlop_params_flint_t fparams;
    fq_default_mat_init_polyvec(tA1, ta1, ctx, mod_ctx);
    fq_default_mat_init_polyvec(tA2, ta2, ctx, mod_ctx);
    fq_default_mat_init_polyvec(tB, tb, ctx, mod_ctx);
    fq_default_mat_init_polyvec(s1, ms1, ctx, mod_ctx);
    fq_default_mat_init_polyvec(m, mm, ctx, mod_ctx);
    fq_default_mat_init_polyvec(s2, ms2, ctx, mod_ctx);
    fq_default_mat_init_polymat(A1, a1, ctx, mod_ctx);
    fq_default_mat_init_polymat(A2prime, a2prime, ctx, mod_ctx);
    fq_default_mat_init_polymat(Bprime, bprime, ctx, mod_ctx);
    abdlop_params_to_flint(fparams, params, ctx, mod_ctx);
    abdlop_commit_flint(ta1, ta2, tb, ms1, mm, ms2, a1, a2prime, bprime, fparams);
}

int abdlop_verify_flint(
    uint8_t hash[32],
    fmpz_mod_poly_t c,
    fq_default_mat_t z1,
    fq_default_mat_t z21,
    fq_default_mat_t h,
    fq_default_mat_t tA1,
    fq_default_mat_t A1,
    fq_default_mat_t A2prime,
    const abdlop_params_flint_t params
) {
    slong kmsis = params->kmsis;
    slong m1 = params->m1;
    slong log2m = params->dcompress->log2m;

    fmpz_mod_poly_t modulus_poly;
    fmpz_mod_poly_init(modulus_poly, params->mod_ctx);
    fq_default_ctx_modulus(modulus_poly, params->ring);
    slong d = fmpz_mod_poly_length(modulus_poly, params->mod_ctx) - 1;
    fmpz_mod_poly_clear(modulus_poly, params->mod_ctx);

    fq_default_mat_t w1, tmp1;
    fq_default_mat_init(w1, kmsis, 1, params->ring);
    fq_default_mat_init(tmp1, kmsis, 1, params->ring);

    fmpz_t l2sqr, bnd, tmp;
    fmpz_init(l2sqr);
    fmpz_init(bnd);
    fmpz_init(tmp);

    fq_default_t cd;
    fq_default_init2(cd, params->ring);

    fq_default_mat_t tmp_mul;
    fq_default_mat_init(tmp_mul, kmsis, 1, params->ring);

    coder_state_t cstate;
    unsigned int outlen;
    uint8_t out[CEIL(log2m * d * kmsis, 8) + 1];
    uint8_t cseed[32];
    fmpz_mod_poly_t c2;
    fmpz_mod_poly_init(c2, params->mod_ctx);

    shake128_state_t hstate;
    int b, accept = 0;

    /* recover w1 */
    fq_default_mat_mul(tmp_mul, A1, z1, params->ring);
    fq_default_mat_set(tmp1, tmp_mul, params->ring);
    fq_default_mat_mul(tmp_mul, A2prime, z21, params->ring);
    fq_default_mat_add(tmp1, tmp1, tmp_mul, params->ring);

    fmpz_mod_poly_t c_poly;
    fmpz_mod_poly_init(c_poly, params->mod_ctx);
    fq_default_get_fmpz_mod_poly(c_poly, cd, params->ring);
    fq_default_set_fmpz_mod_poly(cd, c, params->ring);

    fq_default_mat_scalar_mul(tmp_mul, tA1, cd, params->ring);
    fq_default_mat_sub(tmp1, tmp1, tmp_mul, params->ring);

    fq_default_mat_dcompress_use_ghint(w1, h, tmp1, params);

    /* recover challenge from w1 */
    coder_enc_begin(cstate, out);
    coder_enc_urandom3_flint(
        cstate,
        w1,
        params->ring,
        params->mod_ctx,
        params->dcompress->m,
        log2m
    );
    coder_enc_end(cstate);

    outlen = coder_get_offset(cstate);
    ASSERT_ERR(outlen % 8 == 0);
    ASSERT_ERR(outlen / 8 <= CEIL(log2m * d * kmsis, 8) + 1);
    outlen >>= 3;

    shake128_init(hstate);
    shake128_absorb(hstate, hash, 32);
    shake128_absorb(hstate, out, outlen);
    shake128_squeeze(hstate, cseed, 32);

    fmpz_mod_poly_urand_autostable(c2, d, params->mod_ctx, params->omega, params->log2omega, cseed, 0);

    b = fmpz_mod_poly_equal(c, c2, params->mod_ctx);
    if (!b) goto ret;

    /* check bounds */
    fmpz_set_ui(tmp, 400);
    fmpz_set_ui(bnd, 2 * m1 * d * 962);
    fmpz_mul_2exp(bnd, bnd, 2 * params->log2stdev1);
    fmpz_fdiv_q(bnd, bnd, tmp);

    polyvec_l2sqr_flint(l2sqr, z1, params->ring, params->mod_ctx);
    b = fmpz_cmp(l2sqr, bnd) <= 0;
    if (!b) goto ret;

    fmpz_mod_poly_t gamma_poly;
    fmpz_mod_poly_init(gamma_poly, params->mod_ctx);
    fmpz_mod_poly_set_fmpz(gamma_poly, params->dcompress->gamma, params->mod_ctx);

    fq_default_t g;
    fq_default_init2(g, params->ring);
    fq_default_set_fmpz_mod_poly(g, gamma_poly, params->ring);
    fmpz_mod_poly_clear(gamma_poly, params->mod_ctx);

    fq_default_mat_scalar_mul(tmp1, w1, g, params->ring);
    fq_default_clear(g, params->ring);

    polyvec_l2sqr_flint(l2sqr, tmp1, params->ring, params->mod_ctx);
    b = fmpz_cmp(l2sqr, params->Bsqr) <= 0;
    if (!b) goto ret;

    /* 2*linf(h) <= m */
    fq_default_mat_linf_flint(tmp, h, params->ring, params->mod_ctx);
    b = fmpz_cmp(tmp, params->dcompress->m) <= 0;
    if (!b) goto ret;

    /* update fiat-shamir hash */
    memcpy(hash, cseed, 32);
    accept = 1;

ret:
    shake128_clear(hstate);
    fq_default_mat_clear(w1, params->ring);
    fq_default_mat_clear(tmp1, params->ring);
    fq_default_mat_clear(tmp_mul, params->ring);
    fq_default_clear(cd, params->ring);
    fmpz_clear(l2sqr);
    fmpz_clear(bnd);
    fmpz_clear(tmp);
    fmpz_mod_poly_clear(c2, params->mod_ctx);
    fmpz_mod_poly_clear(c_poly, params->mod_ctx);
    return accept;
}
