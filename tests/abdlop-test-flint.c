#include "abdlop-params1.h"
#include "abdlop-params2.h"
#include "abdlop-params3.h"
#include "abdlop-params4.h"
#include "lazer.h"
#include "test.h"
#include "../src/flint/abdlop.h"
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

static void test_abdlop (uint8_t seed[32], const abdlop_params_t params);

int
main (void)
{
  unsigned int i;
  uint8_t seed[32];

  lazer_init();

  for (i = 0; i < 1; i++)
    {
      printf("\nTEST 1\n");
      bytes_urandom (seed, sizeof (seed));
      test_abdlop (seed, params1);
    }
  // for (i = 0; i < 1; i++)
  //   {
  //     printf("\nTEST 2\n");
  //     bytes_urandom (seed, sizeof (seed));
  //     test_abdlop (seed, params2);
  //   }
  // for (i = 0; i < 1; i++)
  //   {
  //     printf("\nTEST 3\n");
  //    bytes_urandom (seed, sizeof (seed));
  //    test_abdlop (seed, params3);
  //   }
  // for (i = 0; i < 1; i++)
  //   {
  //     printf("\nTEST 4\n");
  //     bytes_urandom (seed, sizeof (seed));
  //     test_abdlop (seed, params4);
  //   }

  printf("\nTEST_PASS\n");
  TEST_PASS ();
}

static void
test_abdlop (uint8_t seed[32], const abdlop_params_t params)
{
  uint8_t hashp[32] = { 0 };
  uint8_t hashv[32] = { 0 };
  uint8_t hashcomm[32] = { 0 };
  polyring_srcptr Rq = params->ring;
  INT_T (lo, Rq->q->nlimbs);
  INT_T (hi, Rq->q->nlimbs);
  uint8_t buf[2];
  uint32_t dom;
  unsigned int i;
  int b;
  polymat_t A1, A2prime, Bprime, A1err, A2primeerr;
  polyvec_t z1err, z21err, herr, tA1err, s1, s2, m, tA1, tA2, tB, z1, z21, h;
  poly_t c;
  const unsigned int l = params->l + params->lext;

  poly_alloc (c, Rq);
  polyvec_alloc (z1err, Rq, params->m1);
  polyvec_alloc (z21err, Rq, params->m2 - params->kmsis);
  polyvec_alloc (herr, Rq, params->kmsis);
  polyvec_alloc (tA1err, Rq, params->kmsis);
  polyvec_alloc (s1, Rq, params->m1);
  polyvec_alloc (s2, Rq, params->m2);
  if (l > 0)
    polyvec_alloc (m, Rq, params->l + params->lext);
  polyvec_alloc (tA1, Rq, params->kmsis);
  polyvec_alloc (tA2, Rq, params->kmsis);
  if (l > 0)
    polyvec_alloc (tB, Rq, params->l + params->lext);
  polyvec_alloc (z1, Rq, params->m1);
  polyvec_alloc (z21, Rq, params->m2 - params->kmsis);
  polyvec_alloc (h, Rq, params->kmsis);
  polymat_alloc (A1, Rq, params->kmsis, params->m1);
  polymat_alloc (A2prime, Rq, params->kmsis, params->m2 - params->kmsis);
  if (l > 0)
    polymat_alloc (Bprime, Rq, params->l + params->lext,
                   params->m2 - params->kmsis);
  polymat_alloc (A1err, Rq, params->kmsis, params->m1);
  polymat_alloc (A2primeerr, Rq, params->kmsis, params->m2 - params->kmsis);

  dom = 0;
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s1, lo, hi, seed, dom++);
  polyvec_urandom_bnd (s2, lo, hi, seed, dom++);
  if (l > 0)
    polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

  /* generate public parameters */

  fq_default_mat_t a1, a2prime, bprime;
  fmpz_mod_ctx_t mod_ctx;
  fq_default_ctx_t ctx;
  abdlop_params_flint_t fparams;
  fq_default_mat_init_polymat(A1, a1, ctx, mod_ctx);
  fq_default_mat_init_polymat(A2prime, a2prime, ctx, mod_ctx);
  fq_default_mat_init_polymat(Bprime, bprime, ctx, mod_ctx);
  abdlop_params_to_flint(fparams, params, ctx, mod_ctx);
  clock_t start = clock();
  abdlop_keygen_flint(a1, a2prime, bprime, seed, fparams);
  clock_t end = clock();
  printf("Flint Key Gen: %f ms\n", (double)(end - start) * 1000.0 / CLOCKS_PER_SEC);

  /* generate proof */

  memcpy (hashp, seed, 32);

  fq_default_mat_t ta1, ta2, tb, ms1, mm, ms2;
  fq_default_mat_init_polyvec(tA1, ta1, ctx, mod_ctx);
  fq_default_mat_init_polyvec(tA2, ta2, ctx, mod_ctx);
  fq_default_mat_init_polyvec(tB, tb, ctx, mod_ctx);
  fq_default_mat_init_polyvec(s1, ms1, ctx, mod_ctx);
  fq_default_mat_init_polyvec(m, mm, ctx, mod_ctx);
  fq_default_mat_init_polyvec(s2, ms2, ctx, mod_ctx);
  abdlop_params_to_flint(fparams, params, ctx, mod_ctx);
  start = clock();
  abdlop_commit_flint(ta1, ta2, tb, ms1, mm, ms2, a1, a2prime, bprime, fparams);
  end = clock();
  printf("Flint Commit: %f ms\n", (double)(end - start) * 1000.0 / CLOCKS_PER_SEC);
  abdlop_hashcomm_flint(hashp, ta1, tb, fparams);

  fmpz_mod_poly_t mc;
  fmpz_mod_poly_init(mc, mod_ctx);
  unsigned int nelems = intvec_get_nelems(poly_get_coeffvec(c));
  for (unsigned int k = 0; k < nelems; k++) {
    fmpz_t cc;
    fmpz_init_int(cc, poly_get_coeff(c, k));
    fmpz_mod_poly_set_coeff_fmpz(mc, k, cc, mod_ctx);
  }

  fq_default_mat_t mz1, mz21, mh;
  fq_default_mat_init_polyvec(z1, mz1, ctx, mod_ctx);
  fq_default_mat_init_polyvec(z21, mz21, ctx, mod_ctx);
  fq_default_mat_init_polyvec(h, mh, ctx, mod_ctx);
  start = clock();
  abdlop_prove_flint(hashp, mc, mz1, mz21, mh, ta2, ms1, ms2, a1, a2prime, seed, fparams);
  end = clock();
  printf("Flint Prove: %f ms\n", (double)(end - start) * 1000.0 / CLOCKS_PER_SEC);

  /* expect successful verification */

  memcpy (hashv, seed, 32);
  abdlop_hashcomm_flint(hashv, ta1, tb, fparams);
  // start = clock();
  // b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
  // end = clock();
  // printf("Lazer Verify: %f ms\n", (double)(end - start) * 1000.0 / CLOCKS_PER_SEC);
  // TEST_EXPECT (b == 1);
  // TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
  //
  // memcpy (hashcomm, seed, 32);
  // abdlop_hashcomm (hashcomm, tA1, tB, params);

  // for (i = 0; i < 10; i++)
  //   {
  //     /* expect verification failures */
  //
  //     bytes_urandom (buf, sizeof (buf));
  //     memcpy (hashv, hashcomm, 32);
  //     hashv[buf[0] % 32] ^= (1 << (buf[1] % 8));
  //     b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
  //     TEST_EXPECT (b == 0);
  //
  //     polyvec_brandom (z1err, 1, seed, dom++);
  //     polyvec_add (z1err, z1err, z1, 0);
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1err, z21, h, tA1, A1, A2prime, params);
  //     TEST_EXPECT (b == 0);
  //
  //     polyvec_brandom (z21err, 1, seed, dom++);
  //     polyvec_add (z21err, z21err, z21, 0);
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1, z21err, h, tA1, A1, A2prime, params);
  //     TEST_EXPECT (b == 0);
  //
  //     polyvec_brandom (herr, 1, seed, dom++);
  //     polyvec_add (herr, herr, h, 0);
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1, z21, herr, tA1, A1, A2prime, params);
  //     TEST_EXPECT (b == 0);
  //
  //     polyvec_brandom (tA1err, 1, seed, dom++); /* sometimes fails XXX */
  //     polyvec_add (tA1err, tA1err, tA1, 0);
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1, z21, h, tA1err, A1, A2prime, params);
  //     TEST_EXPECT (b == 0);
  //
  //     polymat_brandom (A1err, 1, seed, dom++);
  //     polymat_add (A1err, A1err, A1, 0);
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1err, A2prime, params);
  //     TEST_EXPECT (b == 0);
  //
  //     polymat_brandom (A2primeerr, 1, seed, dom++);
  //     polymat_add (A2primeerr, A2primeerr, A2prime, 0);
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2primeerr, params);
  //     TEST_EXPECT (b == 0);
  //
  //     /* expect successful verification */
  //
  //     memcpy (hashv, hashcomm, 32);
  //     b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
  //     TEST_EXPECT (b == 1);
  //     TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
  //   }

  poly_free (c);
  polyvec_free (z1err);
  polyvec_free (z21err);
  polyvec_free (herr);
  polyvec_free (tA1err);
  polyvec_free (s1);
  polyvec_free (s2);
  if (l > 0)
    polyvec_free (m);
  polyvec_free (tA1);
  polyvec_free (tA2);
  if (l > 0)
    polyvec_free (tB);
  polyvec_free (z1);
  polyvec_free (z21);
  polyvec_free (h);
  polymat_free (A1);
  polymat_free (A2prime);
  if (l > 0)
    polymat_free (Bprime);
  polymat_free (A1err);
  polymat_free (A2primeerr);
}
