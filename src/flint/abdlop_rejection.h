#pragma once

#include "abdlop.h"

int rej_standard_flint(
    rng_state_t state,
    const fq_default_mat_t z,
    const fq_default_mat_t v,
    const fq_default_ctx_t ctx,
    const fmpz_mod_ctx_t mod_ctx,
    const fmpz_t scM,
    const fmpz_t sigma2
);

void polyvec_l2sqr_flint(
    fmpz_t r, 
    const fq_default_mat_t a, 
    const fq_default_ctx_t ctx,
    const fmpz_mod_ctx_t mod_ctx
);
