#pragma once

#include "abdlop.h"

void coder_enc_urandom2_flint(
    coder_state_t state,
    fmpz_mod_poly_t v,
    slong degree,
    fmpz_mod_ctx_t ctx,
    const fmpz_t m,
    unsigned int mbits
);

void coder_enc_urandom3_flint(
    coder_state_t state,
    fq_default_mat_t v,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t m,
    unsigned int mbits
);
