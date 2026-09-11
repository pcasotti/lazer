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
    const fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
) {
    fmpz_t r1i, r0i;
    fmpz_init(r1i);
    fmpz_init(r0i);

    slong len = fmpz_mod_poly_length(r, params->mod_ctx);

    for (slong i = 0; i < len; i++) {
        /* r1 = r mod q, in [0,q) (canonical coefficient of r) */
        fmpz_mod_poly_get_coeff_fmpz(r1i, r, i, params->mod_ctx);

        /* r0 = |r1 mod gamma|, in [0,gamma) */
        fmpz_mod(r0i, r1i, params->dcompress->gamma);

        /* center r0 to (-gamma/2, gamma/2] */
        if (fmpz_sgn(r0i) < 0) {
            if (fmpz_cmpabs(r0i, params->dcompress->gammaby2) <= 0) {
                fmpz_add(r0i, r0i, params->dcompress->gamma);
            }
        } else {
            if (fmpz_cmpabs(r0i, params->dcompress->gammaby2) > 0) {
                fmpz_sub(r0i, r0i, params->dcompress->gamma);
            }
        }

        /* r1 = r' - r0 (exact multiple of gamma) */
        fmpz_sub(r1i, r1i, r0i);

        if (fmpz_equal(r1i, params->dcompress->qminus1)) {
            fmpz_zero(r1i);
            fmpz_sub_ui(r0i, r0i, 1);
        } else {
            fmpz_divexact(r1i, r1i, params->dcompress->gamma);
        }

        fmpz_mod_poly_set_coeff_fmpz(r1, i, r1i, params->mod_ctx);
        fmpz_mod_poly_set_coeff_fmpz(r0, i, r0i, params->mod_ctx);
    }

    fmpz_clear(r1i);
    fmpz_clear(r0i);
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

void fmpz_mod_poly_dcompress_make_ghint(
    fmpz_mod_poly_t ret,
    const fmpz_mod_poly_t z,
    const fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
) {
    fmpz_mod_poly_t r1, r0, v0, rz;
    fmpz_mod_poly_init(r1, params->mod_ctx);
    fmpz_mod_poly_init(r0, params->mod_ctx);
    fmpz_mod_poly_init(v0, params->mod_ctx);
    fmpz_mod_poly_init(rz, params->mod_ctx);

    /* rz = r + z (mod q), decompose both like the original */
    fmpz_mod_poly_add(rz, r, z, params->mod_ctx);
    fmpz_mod_poly_dcompress_decompose(r1, r0, r, params);
    fmpz_mod_poly_dcompress_decompose(rz, v0, rz, params);

    /* h = r1(r+z) - r1(r) (mod q) */
    fmpz_mod_poly_sub(ret, rz, r1, params->mod_ctx);

    fmpz_t ret_c, m_half, q_half;
    fmpz_init(ret_c);
    fmpz_init(m_half);
    fmpz_fdiv_q_2exp(m_half, params->dcompress->m, 1);
    fmpz_init(q_half);
    fmpz_fdiv_q_2exp(q_half, params->dcompress->q, 1);

    slong deg = fmpz_mod_poly_length(ret, params->mod_ctx);
    for (slong i = 0; i < deg; i++) {
        fmpz_mod_poly_get_coeff_fmpz(ret_c, ret, i, params->mod_ctx);

        /* signed h, in (-m,m) */
        if (fmpz_cmp(ret_c, q_half) > 0) fmpz_sub(ret_c, ret_c, params->dcompress->q);

        /* h mod m, in [0,m) */
        fmpz_mod(ret_c, ret_c, params->dcompress->m);

        /* center: odd -> [-(m-1)/2,(m-1)/2], even -> (-m/2,m/2] */
        if (fmpz_cmp(ret_c, m_half) > 0) fmpz_sub(ret_c, ret_c, params->dcompress->m);

        fmpz_mod(ret_c, ret_c, params->dcompress->q);
        fmpz_mod_poly_set_coeff_fmpz(ret, i, ret_c, params->mod_ctx);
    }

    fmpz_clear(ret_c);
    fmpz_clear(m_half);
    fmpz_clear(q_half);
    fmpz_mod_poly_clear(r1, params->mod_ctx);
    fmpz_mod_poly_clear(r0, params->mod_ctx);
    fmpz_mod_poly_clear(v0, params->mod_ctx);
    fmpz_mod_poly_clear(rz, params->mod_ctx);
}

void fmpz_mod_poly_dcompress_use_ghint(
    fmpz_mod_poly_t ret,
    const fmpz_mod_poly_t y,
    fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
) {
    fmpz_mod_poly_t r1, r0;
    fmpz_mod_poly_init(r1, params->mod_ctx);
    fmpz_mod_poly_init(r0, params->mod_ctx);

    fmpz_mod_poly_dcompress_decompose(r1, r0, r, params);

    /* ret = r1 + y */
    fmpz_mod_poly_add(ret, r1, y, params->mod_ctx);

    /* reduce each coefficient modulo m, to [0,m) */
    fmpz_t coeff;
    fmpz_init(coeff);
    slong len = fmpz_mod_poly_length(ret, params->mod_ctx);
    for (slong i = 0; i < len; i++) {
        fmpz_mod_poly_get_coeff_fmpz(coeff, ret, i, params->mod_ctx);
        fmpz_mod(coeff, coeff, params->dcompress->m);
        fmpz_mod_poly_set_coeff_fmpz(ret, i, coeff, params->mod_ctx);
    }
    fmpz_clear(coeff);

    fmpz_mod_poly_clear(r1, params->mod_ctx);
    fmpz_mod_poly_clear(r0, params->mod_ctx);
}

void fq_default_mat_dcompress_make_ghint(
    fq_default_mat_t ret,
    const fq_default_mat_t z,
    const fq_default_mat_t r,
    const abdlop_params_flint_t params
) {
    for (slong i = 0; i < fq_default_mat_nrows(ret, params->ring); i++) {
        fq_default_t rete, ze, re;
        fmpz_mod_poly_t retp, zp, rp;
        fq_default_init2(rete, params->ring);
        fq_default_init2(ze, params->ring);
        fq_default_init2(re, params->ring);
        fmpz_mod_poly_init(retp, params->mod_ctx);
        fmpz_mod_poly_init(zp, params->mod_ctx);
        fmpz_mod_poly_init(rp, params->mod_ctx);

        fq_default_mat_entry(rete, ret, i, 0, params->ring);
        fq_default_mat_entry(ze, z, i, 0, params->ring);
        fq_default_mat_entry(re, r, i, 0, params->ring);

        fq_default_get_fmpz_mod_poly(retp, rete, params->ring);
        fq_default_get_fmpz_mod_poly(zp, ze, params->ring);
        fq_default_get_fmpz_mod_poly(rp, re, params->ring);

        fmpz_mod_poly_dcompress_make_ghint(retp, zp, rp, params);

        fq_default_set_fmpz_mod_poly(rete, retp, params->ring);
        fq_default_mat_entry_set(ret, i, 0, rete, params->ring);

        fq_default_clear(rete, params->ring);
        fq_default_clear(ze, params->ring);
        fq_default_clear(re, params->ring);
        fmpz_mod_poly_clear(retp, params->mod_ctx);
        fmpz_mod_poly_clear(zp, params->mod_ctx);
        fmpz_mod_poly_clear(rp, params->mod_ctx);
    }
}

void fq_default_mat_dcompress_use_ghint(
    fq_default_mat_t ret,
    const fq_default_mat_t y,
    const fq_default_mat_t r,
    const abdlop_params_flint_t params
) {
    for (slong i = 0; i < fq_default_mat_nrows(ret, params->ring); i++) {
        fq_default_t rete, ye, re;
        fmpz_mod_poly_t retp, yp, rp;
        fq_default_init2(rete, params->ring);
        fq_default_init2(ye, params->ring);
        fq_default_init2(re, params->ring);
        fmpz_mod_poly_init(retp, params->mod_ctx);
        fmpz_mod_poly_init(yp, params->mod_ctx);
        fmpz_mod_poly_init(rp, params->mod_ctx);

        fq_default_mat_entry(rete, ret, i, 0, params->ring);
        fq_default_mat_entry(ye, y, i, 0, params->ring);
        fq_default_mat_entry(re, r, i, 0, params->ring);

        fq_default_get_fmpz_mod_poly(retp, rete, params->ring);
        fq_default_get_fmpz_mod_poly(yp, ye, params->ring);
        fq_default_get_fmpz_mod_poly(rp, re, params->ring);

        fmpz_mod_poly_dcompress_use_ghint(retp, yp, rp, params);

        fq_default_set_fmpz_mod_poly(rete, retp, params->ring);
        fq_default_mat_entry_set(ret, i, 0, rete, params->ring);

        fq_default_clear(rete, params->ring);
        fq_default_clear(ye, params->ring);
        fq_default_clear(re, params->ring);
        fmpz_mod_poly_clear(retp, params->mod_ctx);
        fmpz_mod_poly_clear(yp, params->mod_ctx);
        fmpz_mod_poly_clear(rp, params->mod_ctx);
    }
}

