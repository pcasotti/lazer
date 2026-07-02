#include "abdlop.h"
#include "src/dom.h"
#include "src/intvec.h"
#include "src/urandom.h"
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
    clock_t start = clock();
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

    clock_t end = clock();
    printf("Flint: %f\n", (double)(end-start));

    // Final Cleanup
    fq_default_mat_window_clear(s21, ctx);
    fq_default_mat_window_clear(s22, ctx);
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
    slong d = params->dcompress->D;

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
    uint8_t out[CEIL (log2m * d * kmsis, 8) + 1];

    shake128_state_t hstate;

    int rej;

    uint32_t dom = 0;
    while (1) {
        // fq_default_mat_grand(y1, params->log2stdev1, yseed, dom);
        dom++;
        // fq_default_mat_grand(y2, params->log2stdev2, yseed, dom);
        dom++;

        fq_default_mat_init_set(w, y22, params->ring);
        fq_default_mat_submul(w, w, A1, y1, params->ring);
        fq_default_mat_submul(w, w, A2prime, y21, params->ring);
        fq_default_mat_dcompress_decompose(w1, w0, w, params);

        coder_enc_begin (cstate, out);
        coder_enc_urandom3_flint(
            cstate,
            w1,
            params->ring,
            params->mod_ctx,
            params->dcompress->m,
            log2m
        );
        coder_enc_end (cstate);

        outlen = coder_get_offset (cstate);
        ASSERT_ERR (outlen % 8 == 0);
        ASSERT_ERR (outlen / 8 <= CEIL (log2m * d * kmsis, 8) + 1);
        outlen >>= 3; /* nbits to nbytes */

        shake128_init (hstate);
        shake128_absorb (hstate, hash, 32);
        shake128_absorb (hstate, out, outlen);
        shake128_squeeze (hstate, cseed, 32);

        // poly_urandom_autostable (c, params->omega, params->log2omega, cseed, 0);

        fq_default_t cd;
        fq_default_init2(cd, params->ring);
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

        fq_default_mat_t tmp;
        fq_default_mat_init_set(tmp, y22, params->ring);
        fq_default_mat_scalar_mul(y22, tA2, cd, params->ring);
        fq_default_mat_sub(y22, y22, tmp, params->ring);
        fq_default_mat_sub(y22, y22, w0, params->ring);
    }
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
        printf("\n"); fq_default_mat_print_pretty(A1, params->ring); printf("\n");
        printf("\n"); fq_default_mat_print_pretty(A2prime, params->ring); printf("\n");
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
        printf("\n"); fq_default_mat_print_pretty(Bprime, params->ring); printf("\n");
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
