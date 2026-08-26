#pragma once

#include "abdlop.h"

void fmpz_mod_poly_power2round(
    fmpz_mod_poly_t ret,
    fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
);

void fq_default_mat_power2round(
    fq_default_mat_t r,
    fq_default_mat_t a,
    const abdlop_params_flint_t params
);

void fmpz_mod_poly_dcompress_decompose(
    fmpz_mod_poly_t r1,
    fmpz_mod_poly_t r0,
    fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
);

void fq_default_mat_dcompress_decompose(
    fq_default_mat_t r1,
    fq_default_mat_t r0,
    fq_default_mat_t r,
    const abdlop_params_flint_t params
);

void fmpz_mod_poly_dcompress_make_ghint(
    fmpz_mod_poly_t ret,
    const fmpz_mod_poly_t z,
    const fmpz_mod_poly_t r,
    const abdlop_params_flint_t params
);

void fq_default_mat_dcompress_make_ghint(
    fq_default_mat_t ret,
    const fq_default_mat_t z,
    const fq_default_mat_t r,
    const abdlop_params_flint_t params
);

