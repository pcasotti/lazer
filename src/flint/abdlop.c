#include "abdlop.h"
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>

void polymat_to_fq_mat(polymat_t src, fq_nmod_mat_t dst, fq_nmod_ctx_t ctx) {
    polymat_fromcrt(src);

    nmod_poly_t ring;
    nmod_poly_init(ring, int_get_i64(src->ring->q));
    nmod_poly_set_coeff_ui(ring, src->ring->d, 1);
    nmod_poly_set_coeff_ui(ring, 0, 1);

    fq_nmod_ctx_init_modulus(ctx, ring, "x");

    fq_nmod_mat_init(
        dst,
        polymat_get_nrows(src),
        polymat_get_ncols(src),
        ctx
    );

    int64_t nrows = fq_nmod_mat_nrows(dst, ctx);
    int64_t ncols = fq_nmod_mat_ncols(dst, ctx);

    for (int64_t i = 0; i < nrows; ++i) {
        for (int64_t j = 0; j < ncols; ++j) {
            nmod_poly_struct *poly = fq_nmod_mat_entry(dst, i, j);
            poly_ptr s = polymat_get_elem(src, i, j);

            int64_t nelems = intvec_get_nelems(poly_get_coeffvec(s));

            for (int64_t k = 0; k < nelems; ++k) {
                ulong c = int_get_i64(poly_get_coeff(s, k));
                nmod_poly_set_coeff_ui(poly, k, c);
            }
        }
    }
}

void polyvec_to_fq_mat(polyvec_t src, fq_nmod_mat_t dst, fq_nmod_ctx_t ctx) {
    polyvec_fromcrt(src);
    unsigned int i, j;
    _VEC_FOREACH_ELEM(src, i) {
        poly_ptr poly = polyvec_get_elem(src, i);
        intvec_ptr coeffs = _get_coeffvec(poly);
        _VEC_FOREACH_ELEM(coeffs, j) {
            int_ptr coeff = intvec_get_elem(coeffs, j);
            int_redp(coeff, coeff, src->ring->q);
        }
    }

    nmod_poly_t ring;
    nmod_poly_init(ring, int_get_i64(src->ring->q));
    nmod_poly_set_coeff_ui(ring, src->ring->d, 1);
    nmod_poly_set_coeff_ui(ring, 0, 1);

    fq_nmod_ctx_init_modulus(ctx, ring, "x");

    fq_nmod_mat_init(
        dst,
        polyvec_get_nelems(src),
        1,
        ctx
    );

    int64_t nrows = fq_nmod_mat_nrows(dst, ctx);

    for (int64_t i = 0; i < nrows; ++i) {
        nmod_poly_t entry;
        nmod_poly_init(entry, int_get_i64(src->ring->q));

        poly_ptr s = polyvec_get_elem(src, i);
        int64_t nelems = intvec_get_nelems(poly_get_coeffvec(s));

        for (int64_t k = 0; k < nelems; ++k) {
            ulong c = int_get_i64(poly_get_coeff(s, k));
            nmod_poly_set_coeff_ui(entry, k, c);
        }

        fq_nmod_mat_entry_set(dst, i, 0, entry, ctx);
    }
}

void abdlop_params_to_flint(
    abdlop_params_flint_t dst,
    const abdlop_params_t src,
    fq_nmod_ctx_t ring
) {
    dst->ring = ring;

    fmpz_set_si(dst->dcompress->q, int_get_i64(src->dcompress->q));
    fmpz_set_si(dst->dcompress->qminus1, int_get_i64(src->dcompress->qminus1));
    fmpz_set_si(dst->dcompress->m, int_get_i64(src->dcompress->m));
    fmpz_set_si(dst->dcompress->mby2, int_get_i64(src->dcompress->mby2));
    fmpz_set_si(dst->dcompress->gamma, int_get_i64(src->dcompress->gamma));
    fmpz_set_si(dst->dcompress->gammaby2, int_get_i64(src->dcompress->gammaby2));
    fmpz_set_si(dst->dcompress->pow2D, int_get_i64(src->dcompress->pow2D));
    fmpz_set_si(dst->dcompress->pow2Dby2, int_get_i64(src->dcompress->pow2Dby2));
    dst->dcompress->D = src->dcompress->D;
    dst->dcompress->m_odd = src->dcompress->m_odd;
    dst->dcompress->log2m = src->dcompress->log2m;

    dst->m1 = src->m1;
    dst->m2 = src->m2;
    dst->l = src->l;
    dst->lext = src->lext;
    dst->kmsis = src->kmsis;
    fmpz_set_si(dst->Bsqr, int_get_i64(src->Bsqr));
    dst->nu = src->nu;
    dst->omega = src->omega;
    dst->log2omega = src->log2omega;
    dst->eta = src->eta;
    dst->rej1 = src->rej1;
    dst->log2stdev1 = src->log2stdev1;
    fmpz_set_si(dst->scM1, int_get_i64(src->scM1));
    fmpz_set_si(dst->stdev1sqr, int_get_i64(src->stdev1sqr));
    dst->rej2 = src->rej2;
    dst->log2stdev2 = src->log2stdev2;
    fmpz_set_si(dst->scM2, int_get_i64(src->scM2));
    fmpz_set_si(dst->stdev2sqr, int_get_i64(src->stdev2sqr));
}

void fq_nmod_mat_addmul(fq_nmod_mat_t r, fq_nmod_mat_t a, fq_nmod_mat_t b, fq_nmod_ctx_t ctx) {
    // TODO `fq_nmod_mat_submul` should work but doesnt somehow so we do this ugly shit instead
    fq_nmod_mat_t tmp;
    fq_nmod_mat_init(tmp, fq_nmod_mat_nrows(a, ctx), 1, ctx);
    fq_nmod_mat_mul(tmp, a, b, ctx);
    fq_nmod_mat_add(r, r, tmp, ctx);
    fq_nmod_mat_clear(tmp, ctx);
}

// FIXME vibe coded
void fq_nmod_mat_dcompress_power2round(
    fq_nmod_mat_t r, 
    fq_nmod_mat_t a, 
    fq_nmod_ctx_t ctx, 
    dcompress_params_flint_t params
) {
    long nrows = fq_nmod_mat_nrows(a, ctx);
    ulong powerD = 1UL << params->D;

    nmod_t mod = ctx->mod;

    nmod_poly_t poly_a, poly_r;
    nmod_poly_init(poly_a, mod.n);
    nmod_poly_init(poly_r, mod.n);

    for (long i = 0; i < nrows; i++) {
        fq_nmod_get_nmod_poly(poly_a, fq_nmod_mat_entry(a, i, 0), ctx);

        nmod_poly_zero(poly_r);

        long deg = nmod_poly_degree(poly_a);
        for (long k = 0; k <= deg; k++) {
            ulong coeff = nmod_poly_get_coeff_ui(poly_a, k);

            // Recompose: (coeff * 2^D) mod q
            ulong res = nmod_mul(coeff, powerD, mod);

            nmod_poly_set_coeff_ui(poly_r, k, res);
        }

        fq_nmod_set_nmod_poly(fq_nmod_mat_entry(r, i, 0), poly_r, ctx);
    }

    nmod_poly_clear(poly_a);
    nmod_poly_clear(poly_r);
}

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
) {
    const unsigned int D = params->dcompress->D;
    const unsigned int kmsis = params->kmsis;
    const unsigned int m1 = params->m1;
    const unsigned int m2 = params->m2;
    const unsigned int l = params->l;
    dcompress_params_flint_struct *dcomp_param = params->dcompress;
    fq_nmod_mat_t s21, s22, tb_, m_, bprime_;

    fq_nmod_mat_window_init(s21, s2, 0, 0, m2 - kmsis, 1, ctx);
    fq_nmod_mat_window_init(s22, s2, m2 - kmsis, 0, m2, 1, ctx);

    if (m1 > 0) {
        fq_nmod_mat_init_set(ta2, s22, ctx);
        fq_nmod_mat_addmul(ta2, a1, s1, ctx);
        fq_nmod_mat_addmul(ta2, a2prime, s21, ctx);

        // TODO
        // fq_nmod_mat_dcompress_power2round(ta1, ta2, ctx, dcomp_param);
        // polyvec_sublshift (tA2, tA1, D);
    }

    if (l > 0) {
        fq_nmod_mat_window_init(m_, m, 0, 0, l, 0, ctx);
        fq_nmod_mat_window_init(tb_, tb, 0, 0, l, 0, ctx);
        fq_nmod_mat_window_init(bprime_, bprime, 0, 0, l, m2 - kmsis, ctx);

        fq_nmod_mat_set(tb_, m_, ctx);
        fq_nmod_mat_submul(tb_, tb_, bprime_, s21, ctx);
    }
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
    fq_nmod_mat_t ta1, ta2, tb, ms1, mm, ms2, a1, a2prime, bprime;
    fq_nmod_ctx_t ctx;
    abdlop_params_flint_t mparams;
    polyvec_to_fq_mat(tA1, ta1, ctx);
    polyvec_to_fq_mat(tA2, ta2, ctx);
    polyvec_to_fq_mat(tB, tb, ctx);
    polyvec_to_fq_mat(s1, ms1, ctx);
    polyvec_to_fq_mat(m, mm, ctx);
    polyvec_to_fq_mat(s2, ms2, ctx);
    polymat_to_fq_mat(A1, a1, ctx);
    polymat_to_fq_mat(A2prime, a2prime, ctx);
    polymat_to_fq_mat(Bprime, bprime, ctx);
    abdlop_params_to_flint(mparams, params, ctx);
    abdlop_commit_flint(ta1, ta2, tb, ms1, mm, ms2, a1, a2prime, bprime, ctx, mparams);
}
