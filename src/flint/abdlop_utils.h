#pragma once

#include "abdlop.h"

void fmpz_init_int(fmpz_t f, int_srcptr n);

void fq_default_mat_init_polyvec(polyvec_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);

void fq_default_mat_init_polymat(polymat_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);

/* flint -> lazer conversions (dst already allocated / inited) */
void fq_default_mat_to_polyvec(polyvec_t dst, const fq_default_mat_t src, polyring_srcptr ring, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);
void fq_default_mat_to_polymat(polymat_t dst, const fq_default_mat_t src, polyring_srcptr ring, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);
/* lazer -> flint (dst already inited in fq_ctx) */
void polyvec_to_fq_default_mat(polyvec_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);

void abdlop_params_to_flint(
    abdlop_params_flint_t dst,
    const abdlop_params_t src,
    fq_default_ctx_t ring,
    fmpz_mod_ctx_t mod_ctx
);
