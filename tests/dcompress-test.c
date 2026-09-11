#include "abdlop-params1.h"
#include "test.h"
#include "src/flint/abdlop.h"

static void test (abdlop_params_srcptr params);
static void test_flint (abdlop_params_srcptr params);

int
main (void)
{
  lazer_init();

  test (params1);
  test_flint (params1);

  TEST_PASS ();
}

static void
test (abdlop_params_srcptr params)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  const unsigned int deg = params->ring->d;
  const unsigned int nlimbs = params->ring->q->nlimbs;
  INTVEC_T (r, deg, nlimbs);
  INTVEC_T (r0, deg, nlimbs);
  INTVEC_T (r1, deg, nlimbs);
  INTVEC_T (r1prime, deg, nlimbs);
  INTVEC_T (z, deg, nlimbs);
  INTVEC_T (y, deg, nlimbs);
  INTVEC_T (yprime, deg, nlimbs);
  INT_T (qminus1, nlimbs);
  INT_T (negqminus1, nlimbs);
  INT_T (gammaby2, nlimbs);
  INT_T (neggammaby2, nlimbs);
  unsigned int i;

  int_set (qminus1, params->dcompress->qminus1);
  int_neg (negqminus1, qminus1);
  int_set (gammaby2, params->dcompress->gammaby2);
  int_neg (neggammaby2, gammaby2);

  for (i = 0; i < 4000; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      dom = 1;
      intvec_urandom_bnd (r, negqminus1, qminus1, seed, dom);
      dom = 2;
      intvec_urandom_bnd (z, neggammaby2, gammaby2, seed, dom);

      dcompress_make_ghint (y, z, r, params->dcompress);
      dcompress_use_ghint (r1prime, y, r, params->dcompress);

      intvec_add (r, r, z);
      intvec_mod (r, r, params->ring->q);
      dcompress_decompose (r1, r0, r, params->dcompress);

      TEST_EXPECT (intvec_eq (r1prime, r1) == 1);
    }
}

/* flint test: mirror round-trip using fmpz_mod_poly_dcompress_* */
static void
test_flint (abdlop_params_srcptr params)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  const unsigned int deg = params->ring->d;
  fq_default_ctx_t ctx;
  fmpz_mod_ctx_t mod_ctx;
  abdlop_params_flint_t fparams;
  fmpz_t q, halfq;
  fmpz_mod_poly_t r, z, y, r1prime, r1, r0, rsum;
  polyvec_t pv;
  fq_default_mat_t minit;
  unsigned int i;

  /* init flint contexts from the lazer ring */
  polyvec_alloc (pv, params->ring, 1);
  fq_default_mat_init_polyvec (pv, minit, ctx, mod_ctx);
  abdlop_params_to_flint (fparams, params, ctx, mod_ctx);
  fq_default_mat_clear (minit, ctx);
  polyvec_free (pv);

  fmpz_init_int (q, params->ring->q);
  fmpz_init (halfq);
  fmpz_fdiv_q_2exp (halfq, q, 1);

  fmpz_mod_poly_init (r, mod_ctx);
  fmpz_mod_poly_init (z, mod_ctx);
  fmpz_mod_poly_init (y, mod_ctx);
  fmpz_mod_poly_init (r1prime, mod_ctx);
  fmpz_mod_poly_init (r1, mod_ctx);
  fmpz_mod_poly_init (r0, mod_ctx);
  fmpz_mod_poly_init (rsum, mod_ctx);

  for (i = 0; i < 4000; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      dom = 1;

      /* r in full modular range, z small (bounded like [-gamma/2,gamma/2]) */
      fmpz_mod_poly_urand (r, deg, mod_ctx, q, params->ring->log2q, seed,
                           dom++);
      fmpz_mod_poly_urand_autostable (z, deg, mod_ctx,
                                      fmpz_get_si (fparams->dcompress
                                                   ->gammaby2),
                                      fparams->dcompress->log2m, seed, dom++);

      fmpz_mod_poly_dcompress_make_ghint (y, z, r, fparams);
      fmpz_mod_poly_dcompress_use_ghint (r1prime, y, r, fparams);

      fmpz_mod_poly_add (rsum, r, z, mod_ctx);
      fmpz_mod_poly_dcompress_decompose (r1, r0, rsum, fparams);

      TEST_EXPECT (fmpz_mod_poly_equal (r1prime, r1, mod_ctx));
    }

  fmpz_mod_poly_clear (rsum, mod_ctx);
  fmpz_mod_poly_clear (r0, mod_ctx);
  fmpz_mod_poly_clear (r1, mod_ctx);
  fmpz_mod_poly_clear (r1prime, mod_ctx);
  fmpz_mod_poly_clear (y, mod_ctx);
  fmpz_mod_poly_clear (z, mod_ctx);
  fmpz_mod_poly_clear (r, mod_ctx);
  fmpz_clear (halfq);
  fmpz_clear (q);
}
