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

void fmpz_init_int(fmpz_t f, int_srcptr n) {
    unsigned int nlimbs = int_get_nlimbs(n);
    fmpz_init2(f, nlimbs);
    ulong f_limbs[nlimbs];
    for (unsigned int i = 0; i < nlimbs; i++) {
        f_limbs[i] = n->limbs[i];
    }
    fmpz_set_ui_array(f, f_limbs, nlimbs);
    if (int_sgn(n) != 1) fmpz_neg(f, f);
}

void fq_default_mat_init_polyvec(polyvec_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx) {
    polyvec_fromcrt(src);

    fmpz_t q;
    fmpz_init_int(q, src->ring->q);

    fmpz_mod_ctx_init(mod_ctx, q);

    fmpz_mod_poly_t ring;
    fmpz_mod_poly_init(ring, mod_ctx);
    fmpz_mod_poly_set_coeff_ui(ring, src->ring->d, 1, mod_ctx);
    fmpz_mod_poly_set_coeff_ui(ring, 0, 1, mod_ctx);

    fq_default_ctx_init_modulus(fq_ctx, ring, mod_ctx, "x");

    fq_default_mat_init(dst, polyvec_get_nelems(src), 1, fq_ctx);

    slong nrows = fq_default_mat_nrows(dst, fq_ctx);
    for (slong i = 0; i < nrows; ++i) {
        poly_ptr s = polyvec_get_elem(src, i);
        unsigned int nelems = intvec_get_nelems(poly_get_coeffvec(s));

        fmpz_mod_poly_t poly;
        fmpz_mod_poly_init(poly, mod_ctx);

        for (unsigned int k = 0; k < nelems; k++) {
            fmpz_t c;
            fmpz_init_int(c, poly_get_coeff(s, k));
            fmpz_mod_poly_set_coeff_fmpz(poly, k, c, mod_ctx);
        }

        fq_default_t entry;
        fq_default_init2(entry, fq_ctx);
        fq_default_set_fmpz_mod_poly(entry, poly, fq_ctx);
        fq_default_mat_entry_set(dst, i, 0, entry, fq_ctx);
    }
}

void fq_default_mat_init_polymat(polymat_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx) {
    polymat_fromcrt(src);

    fmpz_t q;
    fmpz_init_int(q, src->ring->q);

    fmpz_mod_ctx_init(mod_ctx, q);

    fmpz_mod_poly_t ring;
    fmpz_mod_poly_init(ring, mod_ctx);
    fmpz_mod_poly_set_coeff_ui(ring, src->ring->d, 1, mod_ctx);
    fmpz_mod_poly_set_coeff_ui(ring, 0, 1, mod_ctx);

    fq_default_ctx_init_modulus(fq_ctx, ring, mod_ctx, "x");

    fq_default_mat_init(dst, polymat_get_nrows(src), polymat_get_ncols(src), fq_ctx);

    slong nrows = fq_default_mat_nrows(dst, fq_ctx);
    slong ncols = fq_default_mat_ncols(dst, fq_ctx);
    for (slong i = 0; i < nrows; ++i) {
        for (slong j = 0; j < ncols; ++j) {
            poly_ptr s = polymat_get_elem(src, i, j);
            unsigned int nelems = intvec_get_nelems(poly_get_coeffvec(s));

            fmpz_mod_poly_t poly;
            fmpz_mod_poly_init(poly, mod_ctx);

            for (unsigned int k = 0; k < nelems; k++) {
                fmpz_t c;
                fmpz_init_int(c, poly_get_coeff(s, k));
                fmpz_mod_poly_set_coeff_fmpz(poly, k, c, mod_ctx);
            }

            fq_default_t entry;
            fq_default_init2(entry, fq_ctx);
            fq_default_set_fmpz_mod_poly(entry, poly, fq_ctx);
            fq_default_mat_entry_set(dst, i, j, entry, fq_ctx);
        }
    }
}

void abdlop_params_to_flint(
    abdlop_params_flint_t dst,
    const abdlop_params_t src,
    fq_default_ctx_t ring,
    fmpz_mod_ctx_t mod_ctx
) {
    dst->ring = ring;
    dst->mod_ctx = mod_ctx;

    fmpz_init_int(dst->dcompress->q, src->dcompress->q);
    fmpz_init_int(dst->dcompress->qminus1, src->dcompress->qminus1);
    fmpz_init_int(dst->dcompress->m, src->dcompress->m);
    fmpz_init_int(dst->dcompress->mby2, src->dcompress->mby2);
    fmpz_init_int(dst->dcompress->gamma, src->dcompress->gamma);
    fmpz_init_int(dst->dcompress->gammaby2, src->dcompress->gammaby2);
    fmpz_init_int(dst->dcompress->pow2D, src->dcompress->pow2D);
    fmpz_init_int(dst->dcompress->pow2Dby2, src->dcompress->pow2Dby2);
    dst->dcompress->D = src->dcompress->D;
    dst->dcompress->m_odd = src->dcompress->m_odd;
    dst->dcompress->log2m = src->dcompress->log2m;

    dst->m1 = src->m1;
    dst->m2 = src->m2;
    dst->l = src->l;
    dst->lext = src->lext;
    dst->kmsis = src->kmsis;
    fmpz_init_int(dst->Bsqr, src->Bsqr);
    dst->nu = src->nu;
    dst->omega = src->omega;
    dst->log2omega = src->log2omega;
    dst->eta = src->eta;
    dst->rej1 = src->rej1;
    dst->log2stdev1 = src->log2stdev1;
    fmpz_init_int(dst->scM1, src->scM1);
    fmpz_init_int(dst->stdev1sqr, src->stdev1sqr);
    dst->rej2 = src->rej2;
    dst->log2stdev2 = src->log2stdev2;
    fmpz_init_int(dst->scM2, src->scM2);
    fmpz_init_int(dst->stdev2sqr, src->stdev2sqr);
}

void fmpz_mod_poly_power2round(
    fmpz_mod_poly_t ret,
    fmpz_mod_poly_t r,
    abdlop_params_flint_t params
) {
    slong len = fmpz_mod_poly_length(r, params->mod_ctx);

    fmpz_t coeff, high, low, retc, r0, rc;
    fmpz_init(coeff);
    fmpz_init(high);
    fmpz_init(low);
    fmpz_init(retc);
    fmpz_init(r0);
    fmpz_init(rc);

    for (slong i = 0; i < len; i++) {
        fmpz_mod_poly_get_coeff_fmpz(rc, r, i, params->mod_ctx);

        fmpz_mod(retc, rc, params->dcompress->q);
        fmpz_fdiv_r_2exp(r0, retc, params->dcompress->D);

        if (fmpz_cmpabs(r0, params->dcompress->pow2Dby2) > 0) {
            fmpz_sub(r0, r0, params->dcompress->pow2D);
        }

        fmpz_sub(retc, retc, r0);
        fmpz_fdiv_q_2exp(retc, retc, params->dcompress->D);

        fmpz_mul_2exp(r0, retc, params->dcompress->D);
        fmpz_sub(rc, rc, r0);

        fmpz_mod_poly_set_coeff_fmpz(ret, i, retc, params->mod_ctx);
        fmpz_mod_poly_set_coeff_fmpz(r, i, rc, params->mod_ctx);
    }
}

void fq_default_mat_power2round(
    fq_default_mat_t r,
    fq_default_mat_t a,
    abdlop_params_flint_t params
) {
    for (slong i = 0; i < fq_default_mat_nrows(r, params->ring); i++) {
        fq_default_t re, ae;
        fmpz_mod_poly_t rp, ap;
        fq_default_init2(re, params->ring);
        fq_default_init2(ae, params->ring);
        fmpz_mod_poly_init(rp, params->mod_ctx);
        fmpz_mod_poly_init(ap, params->mod_ctx);
        fq_default_mat_entry(re, r, i, 0, params->ring);
        fq_default_get_fmpz_mod_poly(rp, re, params->ring);
        fq_default_mat_entry(ae, a, i, 0, params->ring);
        fq_default_get_fmpz_mod_poly(ap, ae, params->ring);
        fmpz_mod_poly_power2round(rp, ap, params);
        fq_default_set_fmpz_mod_poly(re, rp, params->ring);
        fq_default_mat_entry_set(r, i, 0, re, params->ring);
        fq_default_set_fmpz_mod_poly(ae, ap, params->ring);
        fq_default_mat_entry_set(a, i, 0, ae, params->ring);
    }
}

void fmpz_mod_poly_dcompress_decompose(
    fmpz_mod_poly_t r1,
    fmpz_mod_poly_t r0,
    fmpz_mod_poly_t r,
    abdlop_params_flint_t params
) {
    slong len = fmpz_mod_poly_length(r, params->mod_ctx);

    fmpz_t rc, r0c, r1c, one;
    fmpz_init(rc);
    fmpz_init(r0c);
    fmpz_init(r1c);
    fmpz_init(one);
    fmpz_set_ui(one, 1);

    fmpz_mod_poly_zero(r1, params->mod_ctx);
    fmpz_mod_poly_zero(r0, params->mod_ctx);

    for (slong i = 0; i < len; i++) {
        // Read input coefficient
        fmpz_mod_poly_get_coeff_fmpz(rc, r, i, params->mod_ctx);

        // Map r into standard modular range [0, q-1]
        fmpz_mod(r1c, rc, params->dcompress->q);

        // Extract initial low bits r0 = r1 mod 2^D
        fmpz_fdiv_r_2exp(r0c, r1c, params->dcompress->D);

        // Center r0 into the symmetric interval (-gamma/2, gamma/2]
        if (fmpz_cmpabs(r0c, params->dcompress->pow2Dby2) > 0) {
            fmpz_sub(r0c, r0c, params->dcompress->pow2D);
        }

        // Compute high bits layout: r1 = r1 - r0
        fmpz_sub(r1c, r1c, r0c);

        // Check the boundary wrap-around condition using params->dcompress->q
        if (fmpz_equal(r1c, params->dcompress->q)) { 
            fmpz_zero(r1c);
            fmpz_sub(r0c, r0c, one);
        } else {
            // Divide by gamma (2^D) to get final high bits component
            fmpz_fdiv_q_2exp(r1c, r1c, params->dcompress->D);
        }

        // Map final signed components back to unsigned [0, q-1]
        fmpz_mod(r0c, r0c, params->dcompress->q);
        fmpz_mod(r1c, r1c, params->dcompress->q);

        // Write directly back to output polynomials
        fmpz_mod_poly_set_coeff_fmpz(r1, i, r1c, params->mod_ctx);
        fmpz_mod_poly_set_coeff_fmpz(r0, i, r0c, params->mod_ctx);
    }

    fmpz_clear(rc);
    fmpz_clear(r0c);
    fmpz_clear(r1c);
    fmpz_clear(one);
}

void fq_default_mat_dcompress_decompose(
    fq_default_mat_t r1,
    fq_default_mat_t r0,
    fq_default_mat_t r,
    abdlop_params_flint_t params
) {
    for (slong i = 0; i < fq_default_mat_nrows(r, params->ring); i++) {
        fq_default_t r1e, r0e, re;
        fmpz_mod_poly_t r1p, r0p, rp;
        fq_default_init2(r1e, params->ring);
        fq_default_init2(r0e, params->ring);
        fq_default_init2(re, params->ring);
        fmpz_mod_poly_init(r1p, params->mod_ctx);
        fmpz_mod_poly_init(r0p, params->mod_ctx);
        fmpz_mod_poly_init(rp, params->mod_ctx);
        fq_default_mat_entry(r1e, r1, i, 0, params->ring);
        fq_default_mat_entry(r0e, r0, i, 0, params->ring);
        fq_default_mat_entry(re, r, i, 0, params->ring);
        fq_default_get_fmpz_mod_poly(r1p, r1e, params->ring);
        fq_default_get_fmpz_mod_poly(r0p, r0e, params->ring);
        fq_default_get_fmpz_mod_poly(rp, re, params->ring);
        fmpz_mod_poly_dcompress_decompose(r1p, r0p, rp, params);
        fq_default_set_fmpz_mod_poly(r1e, r1p, params->ring);
        fq_default_set_fmpz_mod_poly(r0e, r0p, params->ring);
        fq_default_set_fmpz_mod_poly(re, rp, params->ring);
        fq_default_mat_entry_set(r1, i, 0, r1e, params->ring);
        fq_default_mat_entry_set(r0, i, 0, r0e, params->ring);
        fq_default_mat_entry_set(r, i, 0, re, params->ring);
    }
}

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
    const unsigned int D = params->dcompress->D;
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

void fmpz_mod_poly_urand(
    fmpz_mod_poly_t r,
    slong len,
    fmpz_mod_ctx_t ctx,
    const fmpz_t mod,
    unsigned int log2mod,
    const uint8_t seed[32],
    uint64_t dom
) {
    slong nlimbs = fmpz_size(mod);
    slong nelems = len;
    INTVEC_T(coeffvec, nelems, nlimbs);

    slong modlen = fmpz_size(mod);
    INT_T(int_mod, modlen);
    ulong m_limbs[nlimbs];
    fmpz_get_signed_ui_array(m_limbs, modlen, mod);
    for (unsigned int i = 0; i < modlen; ++i) {
        int_mod->limbs[i] = m_limbs[i];
    }

    _intvec_urandom(coeffvec, int_mod, log2mod, seed, dom);

    for (slong i = 0; i < nelems; ++i) {
        fmpz_t c;
        fmpz_init_int(c, intvec_get_elem(coeffvec, i));
        fmpz_mod_poly_set_coeff_fmpz(r, i, c, ctx);
    }
}

void fq_default_mat_urand(
    fq_default_mat_t r,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t mod,
    unsigned int log2mod,
    const uint8_t seed[32],
    uint32_t dom
) {
    union dom _dom = { { 0, dom } };

    slong nrows = fq_default_mat_nrows(r, ctx);
    slong ncols = fq_default_mat_ncols(r, ctx);
    for (slong i = 0; i < nrows; ++i) {
        for (slong j = 0; j < ncols; ++j) {
            fq_default_t entry;
            fq_default_init2(entry, ctx);
            fq_default_mat_entry(entry, r, i, j, ctx);

            fmpz_mod_poly_t poly;
            fmpz_mod_poly_init(poly, mod_ctx);
            fq_default_get_fmpz_mod_poly(poly, entry, ctx);

            fmpz_mod_poly_t ring;
            fmpz_mod_poly_init(ring, mod_ctx);
            fq_default_ctx_modulus(ring, ctx);
            slong len = fmpz_mod_poly_length(ring, mod_ctx) -1;

            _dom.d32[0] = i * ncols + j;
            fmpz_mod_poly_urand(poly, len, mod_ctx, mod, log2mod, seed, _dom.d64);

            fq_default_set_fmpz_mod_poly(entry, poly, ctx);
            fq_default_mat_entry_set(r, i, j, entry, ctx);
        }
    }
}

void fq_default_mat_grand(
    fq_default_mat_t r,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t mod,
    unsigned int log2mod,
    const uint8_t seed[32],
    uint32_t dom
) {
    union dom _dom = { { 0, dom } };

    slong nrows = fq_default_mat_nrows(r, ctx);
    slong ncols = fq_default_mat_ncols(r, ctx);
    for (slong i = 0; i < nrows; ++i) {
        for (slong j = 0; j < ncols; ++j) {
            fq_default_t entry;
            fq_default_init2(entry, ctx);
            fq_default_mat_entry(entry, r, i, j, ctx);

            fmpz_mod_poly_t poly;
            fmpz_mod_poly_init(poly, mod_ctx);
            fq_default_get_fmpz_mod_poly(poly, entry, ctx);

            fmpz_mod_poly_t ring;
            fmpz_mod_poly_init(ring, mod_ctx);
            fq_default_ctx_modulus(ring, ctx);
            slong len = fmpz_mod_poly_length(ring, mod_ctx) -1;

            _dom.d32[0] = i * ncols + j;
            fmpz_mod_poly_urand(poly, len, mod_ctx, mod, log2mod, seed, _dom.d64);

            fq_default_set_fmpz_mod_poly(entry, poly, ctx);
            fq_default_mat_entry_set(r, i, j, entry, ctx);
        }
    }
}

static unsigned int uencode_flint(
    uint8_t **byte,
    unsigned int *bit,
    const fmpz_mod_poly_t v,
    const fmpz_mod_ctx_t ctx,
    const fmpz_t m,
    unsigned int mbits
) {
    unsigned int nbits = 0;
    unsigned int nbytes = 0;

    uint8_t *_byte = *byte;
    unsigned int _bit = *bit;

    _byte[0] &= ~((uint8_t)(~0) << _bit);

    slong len = fmpz_mod_poly_degree(v, ctx);
    fmpz_t elem;
    fmpz_init(elem);

    for (slong i = 0; i < len; i++) {
        fmpz_mod_poly_get_coeff_fmpz(elem, v, i, ctx);

        for (unsigned int k = 0; k < mbits; k++) {
            int bit_val = fmpz_tstbit(elem, k);

            _byte[nbytes] |= ((uint8_t)bit_val << _bit);

            _bit++;
            if (_bit == 8) {
                _bit = 0;
                nbytes++;
                _byte[nbytes] = 0;
            }
            nbits++;
        }
    }

    *byte = _byte + nbytes;
    *bit = _bit;

    fmpz_clear(elem);
    return nbits;
}

void coder_enc_urandom2_flint(
    coder_state_t state,
    fmpz_mod_poly_t v,
    fmpz_mod_ctx_t ctx,
    const fmpz_t m,
    unsigned int mbits
) {
    unsigned int nbits = uencode_flint(
        &(state->out),
        &(state->bit_off),
        v,
        ctx,
        m,
        mbits
    );
    state-> byte_off += (nbits >> 3);
}

void coder_enc_urandom3_flint(
    coder_state_t state,
    fq_default_mat_t v,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t m,
    unsigned int mbits
) {
    for (slong i = 0; i < fq_default_mat_nrows(v, ctx); i++) {
        fq_default_t elem;
        fmpz_mod_poly_t poly;
        fq_default_init2(elem, ctx);
        fmpz_mod_poly_init(poly, mod_ctx);
        fq_default_mat_entry(elem, v, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(poly, elem, ctx);
        coder_enc_urandom2_flint(state, poly, mod_ctx, m, mbits);
        fq_default_set_fmpz_mod_poly(elem, poly, ctx);
        fq_default_mat_entry_set(v, i, 0, elem, ctx);
    }
}

static void rng_urandom_fmpz(fmpz_t out, rng_state_t state) {
    uint8_t bytes[16]; // 128 bits
    rng_urandom(state, bytes, 16);

    mpz_t gmp_num;
    mpz_init(gmp_num);

    // GMP native import: 
    // 16 items, MSB first (1), 1 byte size (1), native endianness (0), 0 nails
    mpz_import(gmp_num, 16, 1, 1, 0, 0, bytes);

    // Safely assign the GMP integer to the FLINT fmpz_t
    fmpz_set_mpz(out, gmp_num);

    mpz_clear(gmp_num);
}

int rej_standard_flint(
    rng_state_t state,
    const fq_default_mat_t z,
    const fq_default_mat_t v,
    const fq_default_ctx_t ctx,
    const fmpz_mod_ctx_t mod_ctx,
    const fmpz_t scM,
    const fmpz_t sigma2
) {
    int reject;
    slong nrows = fq_default_mat_nrows(z, ctx);

    fmpz_t sigma2dbl, t1, t2, u, mu, dot_tmp;
    fmpz_init(sigma2dbl);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(u);
    fmpz_init(mu);
    fmpz_init(dot_tmp);

    // Fix: Correct initialization using fmpz_set on the modulus pointer
    fmpz_t q, q_half;
    fmpz_init(q);
    fmpz_init(q_half);
    fmpz_set(q, fmpz_mod_ctx_modulus(mod_ctx));
    fmpz_fdiv_q_2exp(q_half, q, 1); // q / 2

    // Sample 128-bit random token using project's native RNG layout
    rng_urandom_fmpz(u, state);
    fmpz_mul(mu, scM, u);

    fmpz_zero(t1);
    fmpz_zero(t2);
    
    fmpz_mod_poly_t p_z, p_v;
    fmpz_mod_poly_init(p_z, mod_ctx);
    fmpz_mod_poly_init(p_v, mod_ctx);

    fq_default_t entry_z, entry_v;
    fq_default_init2(entry_z, ctx);
    fq_default_init2(entry_v, ctx);

    // Compute standard Euclidean geometric dot products over the vector coordinates
    for (slong i = 0; i < nrows; i++) {
        fq_default_mat_entry(entry_z, z, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(p_z, entry_z, ctx);

        fq_default_mat_entry(entry_v, v, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(p_v, entry_v, ctx);

        slong len_z = fmpz_mod_poly_length(p_z, mod_ctx);
        for (slong k = 0; k < len_z; k++) {
            fmpz_t c_z, c_v;
            fmpz_init(c_z); fmpz_init(c_v);
            fmpz_mod_poly_get_coeff_fmpz(c_z, p_z, k, mod_ctx);
            fmpz_mod_poly_get_coeff_fmpz(c_v, p_v, k, mod_ctx);

            // Center coefficients from unsigned [0, q-1] into signed [-(q-1)/2, (q-1)/2]
            if (fmpz_cmp(c_z, q_half) > 0) fmpz_sub(c_z, c_z, q);
            if (fmpz_cmp(c_v, q_half) > 0) fmpz_sub(c_v, c_v, q);

            fmpz_mul(dot_tmp, c_z, c_v);
            fmpz_add(t1, t1, dot_tmp);
            fmpz_clear(c_z); fmpz_clear(c_v);
        }

        slong len_v = fmpz_mod_poly_length(p_v, mod_ctx);
        for (slong k = 0; k < len_v; k++) {
            fmpz_t c_v;
            fmpz_init(c_v);
            fmpz_mod_poly_get_coeff_fmpz(c_v, p_v, k, mod_ctx);

            // Center coefficients from unsigned [0, q-1] into signed [-(q-1)/2, (q-1)/2]
            if (fmpz_cmp(c_v, q_half) > 0) fmpz_sub(c_v, c_v, q);

            fmpz_mul(dot_tmp, c_v, c_v);
            fmpz_add(t2, t2, dot_tmp);
            fmpz_clear(c_v);
        }
    }

    // Set up standard exponential targets: t2 = <v,v> - 2*<z,v>
    fmpz_mul_2exp(t1, t1, 1);       
    fmpz_sub(t2, t2, t1);           
    fmpz_mul_2exp(sigma2dbl, sigma2, 1); 

    // --- MPFR Core Evaluation Block (Exact Lazer Equivalent via GMP casting) ---
    mpfr_t nom, denom, t3, t4;
    mpfr_init2(nom, 128);
    mpfr_init2(denom, 128);
    mpfr_init2(t3, 128);
    mpfr_init2(t4, 128);

    // Explicitly unwrap fmpz_t to mpz_t to bypass fragile conditionally-compiled macros
    mpz_t gmp_nom, gmp_denom, gmp_mu;
    mpz_init(gmp_nom);
    mpz_init(gmp_denom);
    mpz_init(gmp_mu);

    fmpz_get_mpz(gmp_nom, t2);
    fmpz_get_mpz(gmp_denom, sigma2dbl);
    fmpz_get_mpz(gmp_mu, mu);

    // Cast the mpz integers safely into 128-bit float structures using round-to-nearest
    mpfr_set_z(nom, gmp_nom, MPFR_RNDN);
    mpfr_set_z(denom, gmp_denom, MPFR_RNDN);
    mpfr_set_z(t4, gmp_mu, MPFR_RNDN);

    // Compute the target distribution exponent threshold
    mpfr_div(t3, nom, denom, MPFR_RNDN);   // (-2<z,v> + <v,v>) / (2*sigma^2)
    mpfr_exp(t3, t3, MPFR_RNDN);           // exp(...)
    mpfr_mul_2exp(t3, t3, 256, MPFR_RNDN); // 2^256 * exp(...)

    // Perform the strict validation evaluation boundary comparison
    if (mpfr_cmp(t4, t3) > 0) {
        reject = 1;
    } else {
        reject = 0;
    }

    // Clean up local MPFR and intermediate GMP allocations
    mpz_clear(gmp_nom); mpz_clear(gmp_denom); mpz_clear(gmp_mu);
    mpfr_clear(nom); mpfr_clear(denom); mpfr_clear(t3); mpfr_clear(t4);

    // Clean up working fmpz allocations and context registers
    fmpz_clear(q); fmpz_clear(q_half);
    fmpz_clear(sigma2dbl); fmpz_clear(t1); fmpz_clear(t2);
    fmpz_clear(u); fmpz_clear(mu); fmpz_clear(dot_tmp);
    fmpz_mod_poly_clear(p_z, mod_ctx);
    fmpz_mod_poly_clear(p_v, mod_ctx);
    fq_default_clear(entry_z, ctx);
    fq_default_clear(entry_v, ctx);

    return reject;
}

void polyvec_l2sqr_flint(
    fmpz_t r, 
    const fq_default_mat_t a, 
    const fq_default_ctx_t ctx,
    const fmpz_mod_ctx_t mod_ctx
) {
    slong nrows = fq_default_mat_nrows(a, ctx);
    
    fmpz_set_ui(r, 0); // Initialize accumulator r = 0
    
    // Set up modulus details to center coefficients correctly around zero
    fmpz_t q, q_half, c_sqr, centered_c;
    fmpz_init(q);
    fmpz_init(q_half);
    fmpz_init(c_sqr);
    fmpz_init(centered_c);
    
    fmpz_set(q, fmpz_mod_ctx_modulus(mod_ctx));
    fmpz_fdiv_q_2exp(q_half, q, 1); // q_half = q / 2

    // Temporary workspace structures
    fmpz_mod_poly_t poly_elem;
    fmpz_mod_poly_init(poly_elem, mod_ctx);
    
    fq_default_t mat_entry;
    fq_default_init2(mat_entry, ctx);

    // Loop through each polynomial element in the vector (column 0)
    for (slong i = 0; i < nrows; i++) 
    {
        // 1. Extract the polynomial from the matrix row
        fq_default_mat_entry(mat_entry, a, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(poly_elem, mat_entry, ctx);
        
        // 2. Loop through every coefficient inside the polynomial
        slong deg = fmpz_mod_poly_length(poly_elem, mod_ctx);
        for (slong j = 0; j < deg; j++) 
        {
            fmpz_mod_poly_get_coeff_fmpz(centered_c, poly_elem, j, mod_ctx);
            
            // 3. Center the unsigned coefficient to the signed range: [-(q-1)/2, (q-1)/2]
            if (fmpz_cmp(centered_c, q_half) > 0) 
            {
                fmpz_sub(centered_c, centered_c, q);
            }
            
            // 4. Square the centered coefficient: c_sqr = centered_c^2
            fmpz_mul(c_sqr, centered_c, centered_c);
            
            // 5. Accumulate into total norm: r = r + c_sqr
            fmpz_add(r, r, c_sqr);
        }
    }

    // Memory cleanup
    fmpz_clear(q);
    fmpz_clear(q_half);
    fmpz_clear(c_sqr);
    fmpz_clear(centered_c);
    fmpz_mod_poly_clear(poly_elem, mod_ctx);
    fq_default_clear(mat_entry, ctx);
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
    abdlop_params_t params
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
