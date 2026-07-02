#include "abdlop_coder.h"
#include <flint/flint.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>

static unsigned int uencode_flint(
    uint8_t **byte,
    unsigned int *bit,
    const fmpz_mod_poly_t v,
    const fmpz_mod_ctx_t ctx,
    const fmpz_t m,
    unsigned int mbits
) {
    unsigned int nbits = 0;
    unsigned int nbytes = 0;

    uint8_t *_byte = *byte;
    unsigned int _bit = *bit;

    _byte[0] &= ~((uint8_t)(~0) << _bit);

    slong len = fmpz_mod_poly_degree(v, ctx);
    fmpz_t elem;
    fmpz_init(elem);

    for (slong i = 0; i < len; i++) {
        fmpz_mod_poly_get_coeff_fmpz(elem, v, i, ctx);

        for (unsigned int k = 0; k < mbits; k++) {
            int bit_val = fmpz_tstbit(elem, k);

            _byte[nbytes] |= ((uint8_t)bit_val << _bit);

            _bit++;
            if (_bit == 8) {
                _bit = 0;
                nbytes++;
                _byte[nbytes] = 0;
            }
            nbits++;
        }
    }

    *byte = _byte + nbytes;
    *bit = _bit;

    fmpz_clear(elem);
    return nbits;
}

void coder_enc_urandom2_flint(
    coder_state_t state,
    fmpz_mod_poly_t v,
    fmpz_mod_ctx_t ctx,
    const fmpz_t m,
    unsigned int mbits
) {
    unsigned int nbits = uencode_flint(
        &(state->out),
        &(state->bit_off),
        v,
        ctx,
        m,
        mbits
    );
    state-> byte_off += (nbits >> 3);
}

void coder_enc_urandom3_flint(
    coder_state_t state,
    fq_default_mat_t v,
    fq_default_ctx_t ctx,
    fmpz_mod_ctx_t mod_ctx,
    const fmpz_t m,
    unsigned int mbits
) {
    for (slong i = 0; i < fq_default_mat_nrows(v, ctx); i++) {
        fq_default_t elem;
        fmpz_mod_poly_t poly;
        fq_default_init2(elem, ctx);
        fmpz_mod_poly_init(poly, mod_ctx);
        fq_default_mat_entry(elem, v, i, 0, ctx);
        fq_default_get_fmpz_mod_poly(poly, elem, ctx);
        coder_enc_urandom2_flint(state, poly, mod_ctx, m, mbits);
        fq_default_set_fmpz_mod_poly(elem, poly, ctx);
        fq_default_mat_entry_set(v, i, 0, elem, ctx);
    }
}
