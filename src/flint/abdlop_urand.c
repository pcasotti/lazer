#include "abdlop_urand.h"
#include "src/dom.h"
#include "src/intvec.h"
#include "src/urandom.h"
#include <flint/flint.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>

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

void fmpz_mod_poly_grand(
    fmpz_mod_poly_t r,
    slong len,
    fmpz_mod_ctx_t ctx,
    unsigned int log2o,
    const uint8_t seed[32],
    uint64_t dom
) {
    INTVEC_T(coeffvec, len, 2); // 2 limbs is generally enough for Gaussian noise, but let's make it robust

    _intvec_grandom(coeffvec, log2o, seed, dom);

    for (slong i = 0; i < len; ++i) {
        fmpz_t c;
        fmpz_init_int(c, intvec_get_elem(coeffvec, i));
        fmpz_mod_poly_set_coeff_fmpz(r, i, c, ctx);
    }
}

void fmpz_mod_poly_urand_autostable(
    fmpz_mod_poly_t r,
    slong len,
    fmpz_mod_ctx_t ctx,
    int64_t bnd,
    unsigned int log2,
    const uint8_t seed[32],
    uint64_t dom
) {
    const unsigned int nchal_coeffs = len / 2;
    unsigned int i;
    int64_t chal_coeffs[nchal_coeffs];
    int64_t mod = 2 * bnd + 1;

    ASSERT_ERR (len % 2 == 0);

    /* coeffs in [0, 2*bnd] */
    _urandom_i64 (chal_coeffs, nchal_coeffs, mod, log2, seed, dom);

    fmpz_t c;
    fmpz_init(c);

    for (i = 0; i < nchal_coeffs; i++)
    {
        /* coeffs in [-bnd, bnd] */
        int64_t val = chal_coeffs[i] - bnd;
        fmpz_set_si(c, val);
        fmpz_mod_poly_set_coeff_fmpz(r, i, c, ctx);
    }
    fmpz_zero(c);
    fmpz_mod_poly_set_coeff_fmpz(r, nchal_coeffs, c, ctx);
    for (i = nchal_coeffs + 1; i < len; i++)
    {
        int64_t val = -chal_coeffs[2 * nchal_coeffs - i];
        fmpz_set_si(c, val);
        fmpz_mod_poly_set_coeff_fmpz(r, i, c, ctx);
    }
    fmpz_clear(c);
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
    (void)mod;
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
            slong len = fmpz_mod_poly_length(ring, mod_ctx) - 1;

            _dom.d32[0] = i * ncols + j;
            fmpz_mod_poly_grand(poly, len, mod_ctx, log2mod, seed, _dom.d64);

            fq_default_set_fmpz_mod_poly(entry, poly, ctx);
            fq_default_mat_entry_set(r, i, j, entry, ctx);
        }
    }
}
