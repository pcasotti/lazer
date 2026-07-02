#include "abdlop_utils.h"
#include "src/dom.h"
#include "src/intvec.h"
#include "src/urandom.h"
#include <flint/flint.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_types.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod.h>

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
