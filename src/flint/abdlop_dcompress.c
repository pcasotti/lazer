#include "abdlop_dcompress.h"
#include <flint/flint.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod.h>

void fmpz_mod_poly_power2round(
    fmpz_mod_poly_t ret,
    fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
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
    const abdlop_params_flint_t params
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
    const abdlop_params_flint_t params
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
    const abdlop_params_flint_t params
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
