#include "abdlop_rejection.h"
#include <mpfr.h>
#include <flint/flint.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod.h>

static void rng_urandom_fmpz(fmpz_t out, rng_state_t state) {
    uint8_t bytes[16]; // 128 bits
    rng_urandom(state, bytes, 16);

    mpz_t gmp_num;
    mpz_init(gmp_num);

    // GMP native import: 
    // 16 items, MSB first (1), 1 byte size (1), native endianness (0), 0 nails
    mpz_import(gmp_num, 16, 1, 1, 0, 0, bytes);

    // Safely assign the GMP integer to the FLINT fmpz_t
    fmpz_set_mpz(out, gmp_num);

    mpz_clear(gmp_num);
}

int rej_standard_flint(
    rng_state_t state,
    const fq_default_mat_t z,
    const fq_default_mat_t v,
    const fq_default_ctx_t ctx,
    const fmpz_mod_ctx_t mod_ctx,
    const fmpz_t scM,
    const fmpz_t sigma2
) {
    int reject;
    slong nrows = fq_default_mat_nrows(z, ctx);

    fmpz_t sigma2dbl, t1, t2, u, mu, dot_tmp;
    fmpz_init(sigma2dbl);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(u);
    fmpz_init(mu);
    fmpz_init(dot_tmp);

    // Fix: Correct initialization using fmpz_set on the modulus pointer
    fmpz_t q, q_half;
    fmpz_init(q);
    fmpz_init(q_half);
    fmpz_set(q, fmpz_mod_ctx_modulus(mod_ctx));
    fmpz_fdiv_q_2exp(q_half, q, 1); // q / 2

    // Sample 128-bit random token using project's native RNG layout
    rng_urandom_fmpz(u, state);
    fmpz_mul(mu, scM, u);

    fmpz_zero(t1);
    fmpz_zero(t2);
    
    fmpz_mod_poly_t p_z, p_v;
    fmpz_mod_poly_init(p_z, mod_ctx);
    fmpz_mod_poly_init(p_v, mod_ctx);

    fq_default_t entry_z, entry_v;
    fq_default_init2(entry_z, ctx);
    fq_default_init2(entry_v, ctx);

    // Compute standard Euclidean geometric dot products over the vector coordinates
    for (slong i = 0; i < nrows; i++) {
        fq_default_mat_entry(entry_z, z, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(p_z, entry_z, ctx);

        fq_default_mat_entry(entry_v, v, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(p_v, entry_v, ctx);

        slong len_z = fmpz_mod_poly_length(p_z, mod_ctx);
        for (slong k = 0; k < len_z; k++) {
            fmpz_t c_z, c_v;
            fmpz_init(c_z); fmpz_init(c_v);
            fmpz_mod_poly_get_coeff_fmpz(c_z, p_z, k, mod_ctx);
            fmpz_mod_poly_get_coeff_fmpz(c_v, p_v, k, mod_ctx);

            // Center coefficients from unsigned [0, q-1] into signed [-(q-1)/2, (q-1)/2]
            if (fmpz_cmp(c_z, q_half) > 0) fmpz_sub(c_z, c_z, q);
            if (fmpz_cmp(c_v, q_half) > 0) fmpz_sub(c_v, c_v, q);

            fmpz_mul(dot_tmp, c_z, c_v);
            fmpz_add(t1, t1, dot_tmp);
            fmpz_clear(c_z); fmpz_clear(c_v);
        }

        slong len_v = fmpz_mod_poly_length(p_v, mod_ctx);
        for (slong k = 0; k < len_v; k++) {
            fmpz_t c_v;
            fmpz_init(c_v);
            fmpz_mod_poly_get_coeff_fmpz(c_v, p_v, k, mod_ctx);

            // Center coefficients from unsigned [0, q-1] into signed [-(q-1)/2, (q-1)/2]
            if (fmpz_cmp(c_v, q_half) > 0) fmpz_sub(c_v, c_v, q);

            fmpz_mul(dot_tmp, c_v, c_v);
            fmpz_add(t2, t2, dot_tmp);
            fmpz_clear(c_v);
        }
    }

    // Set up standard exponential targets: t2 = <v,v> - 2*<z,v>
    fmpz_mul_2exp(t1, t1, 1);       
    fmpz_sub(t2, t2, t1);           
    fmpz_mul_2exp(sigma2dbl, sigma2, 1); 

    // --- MPFR Core Evaluation Block (Exact Lazer Equivalent via GMP casting) ---
    mpfr_t nom, denom, t3, t4;
    mpfr_init2(nom, 128);
    mpfr_init2(denom, 128);
    mpfr_init2(t3, 128);
    mpfr_init2(t4, 128);

    // Explicitly unwrap fmpz_t to mpz_t to bypass fragile conditionally-compiled macros
    mpz_t gmp_nom, gmp_denom, gmp_mu;
    mpz_init(gmp_nom);
    mpz_init(gmp_denom);
    mpz_init(gmp_mu);

    fmpz_get_mpz(gmp_nom, t2);
    fmpz_get_mpz(gmp_denom, sigma2dbl);
    fmpz_get_mpz(gmp_mu, mu);

    // Cast the mpz integers safely into 128-bit float structures using round-to-nearest
    mpfr_set_z(nom, gmp_nom, MPFR_RNDN);
    mpfr_set_z(denom, gmp_denom, MPFR_RNDN);
    mpfr_set_z(t4, gmp_mu, MPFR_RNDN);

    // Compute the target distribution exponent threshold
    mpfr_div(t3, nom, denom, MPFR_RNDN);   // (-2<z,v> + <v,v>) / (2*sigma^2)
    mpfr_exp(t3, t3, MPFR_RNDN);           // exp(...)
    mpfr_mul_2exp(t3, t3, 256, MPFR_RNDN); // 2^256 * exp(...)

    // Perform the strict validation evaluation boundary comparison
    if (mpfr_cmp(t4, t3) > 0) {
        reject = 1;
    } else {
        reject = 0;
    }

    // Clean up local MPFR and intermediate GMP allocations
    mpz_clear(gmp_nom); mpz_clear(gmp_denom); mpz_clear(gmp_mu);
    mpfr_clear(nom); mpfr_clear(denom); mpfr_clear(t3); mpfr_clear(t4);

    // Clean up working fmpz allocations and context registers
    fmpz_clear(q); fmpz_clear(q_half);
    fmpz_clear(sigma2dbl); fmpz_clear(t1); fmpz_clear(t2);
    fmpz_clear(u); fmpz_clear(mu); fmpz_clear(dot_tmp);
    fmpz_mod_poly_clear(p_z, mod_ctx);
    fmpz_mod_poly_clear(p_v, mod_ctx);
    fq_default_clear(entry_z, ctx);
    fq_default_clear(entry_v, ctx);

    return reject;
}

void polyvec_l2sqr_flint(
    fmpz_t r, 
    const fq_default_mat_t a, 
    const fq_default_ctx_t ctx,
    const fmpz_mod_ctx_t mod_ctx
) {
    slong nrows = fq_default_mat_nrows(a, ctx);
    
    fmpz_set_ui(r, 0); // Initialize accumulator r = 0
    
    // Set up modulus details to center coefficients correctly around zero
    fmpz_t q, q_half, c_sqr, centered_c;
    fmpz_init(q);
    fmpz_init(q_half);
    fmpz_init(c_sqr);
    fmpz_init(centered_c);
    
    fmpz_set(q, fmpz_mod_ctx_modulus(mod_ctx));
    fmpz_fdiv_q_2exp(q_half, q, 1); // q_half = q / 2

    // Temporary workspace structures
    fmpz_mod_poly_t poly_elem;
    fmpz_mod_poly_init(poly_elem, mod_ctx);
    
    fq_default_t mat_entry;
    fq_default_init2(mat_entry, ctx);

    // Loop through each polynomial element in the vector (column 0)
    for (slong i = 0; i < nrows; i++) 
    {
        // 1. Extract the polynomial from the matrix row
        fq_default_mat_entry(mat_entry, a, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(poly_elem, mat_entry, ctx);
        
        // 2. Loop through every coefficient inside the polynomial
        slong deg = fmpz_mod_poly_length(poly_elem, mod_ctx);
        for (slong j = 0; j < deg; j++) 
        {
            fmpz_mod_poly_get_coeff_fmpz(centered_c, poly_elem, j, mod_ctx);
            
            // 3. Center the unsigned coefficient to the signed range: [-(q-1)/2, (q-1)/2]
            if (fmpz_cmp(centered_c, q_half) > 0) 
            {
                fmpz_sub(centered_c, centered_c, q);
            }
            
            // 4. Square the centered coefficient: c_sqr = centered_c^2
            fmpz_mul(c_sqr, centered_c, centered_c);
            
            // 5. Accumulate into total norm: r = r + c_sqr
            fmpz_add(r, r, c_sqr);
        }
    }

    // Memory cleanup
    fmpz_clear(q);
    fmpz_clear(q_half);
    fmpz_clear(c_sqr);
    fmpz_clear(centered_c);
    fmpz_mod_poly_clear(poly_elem, mod_ctx);
    fq_default_clear(mat_entry, ctx);
}
