#pragma once

#include "abdlop.h"

void fmpz_mod_poly_urand(
    fmpz_mod_poly_t r,
    slong len,
    fmpz_mod_ctx_t ctx,
    const fmpz_t mod,
    unsigned int log2mod,
    const uint8_t seed[32],
    uint64_t dom
);

void fq_default_mat_urand(
    fq_default_mat_t r,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t mod,
    unsigned int log2mod,
    const uint8_t seed[32],
    uint32_t dom
);

void fq_default_mat_grand(
    fq_default_mat_t r,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t mod,
    unsigned int log2mod,
    const uint8_t seed[32],
    uint32_t dom
);
