#pragma once

#include "abdlop.h"

void fmpz_init_int(fmpz_t f, int_srcptr n);

void fq_default_mat_init_polyvec(polyvec_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);

void fq_default_mat_init_polymat(polymat_t src, fq_default_mat_t dst, fq_default_ctx_t fq_ctx, fmpz_mod_ctx_t mod_ctx);

void abdlop_params_to_flint(
    abdlop_params_flint_t dst,
    const abdlop_params_t src,
    fq_default_ctx_t ring,
    fmpz_mod_ctx_t mod_ctx
);
