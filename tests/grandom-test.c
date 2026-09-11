#include "abdlop-params1.h"
#include "test.h"
#include "src/flint/abdlop.h"
#include <math.h>

#define SP_NSAMPLES 5000000
#define SP_LOG2O 10

#define MP_LOG2NSAMPLES 22
#define MP_LOG2O 128

#define SAMPLEVEC_DIM 4

#define FLINT_NSAMPLES 5000
#define FLINT_LOG2O 10

static void single_prec (void);
static void multiple_prec (void);
static void flint_prec (void);

int
main (void)
{
  lazer_init();

  single_prec ();
  multiple_prec ();
  flint_prec ();
  TEST_PASS ();
}

/* flint test: fmpz_mod_poly_grand distribution (center-reduced coefficients) */
static void
flint_prec (void)
{
  const long double sigma = (long double)1.55 * (1 << FLINT_LOG2O);
  polyring_srcptr Rq = params1->ring;
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  fmpz_t mod;
  fmpz_mod_ctx_t ctx;
  fmpz_mod_poly_t poly;
  fmpz_t q;
  fmpz_t halfq;
  fmpz_t coeff;
  long double sum = 0, sum_sqr = 0;
  long double empiric_mean, empiric_var;
  unsigned int i;

  bytes_urandom (seed, sizeof (seed));

  fmpz_init_int (mod, Rq->q);
  fmpz_init (q);
  fmpz_init (halfq);
  fmpz_init (coeff);
  fmpz_set (q, mod);
  fmpz_fdiv_q_2exp (halfq, q, 1);

  fmpz_mod_ctx_init (ctx, mod);
  fmpz_mod_poly_init (poly, ctx);

  for (i = 0; i < FLINT_NSAMPLES; i++)
    {
      slong k;
      int64_t v;
      fmpz_mod_poly_grand (poly, SAMPLEVEC_DIM, ctx, mod, FLINT_LOG2O, seed,
                           dom);
      dom++;

      for (k = 0; k < SAMPLEVEC_DIM; k++)
        {
          fmpz_mod_poly_get_coeff_fmpz (coeff, poly, k, ctx);
          /* center-reduce into [-q/2, q/2) */
          if (fmpz_cmp (coeff, halfq) > 0)
            fmpz_sub (coeff, coeff, q);
          v = fmpz_get_si (coeff);
          sum += v;
          sum_sqr += (long double)v * v;
        }
    }

  empiric_mean = sum / (FLINT_NSAMPLES * SAMPLEVEC_DIM);
  empiric_var = sum_sqr / (FLINT_NSAMPLES * SAMPLEVEC_DIM);

  printf ("flint gaussian test\n");
  printf ("expected variance: %Lf\n", sigma * sigma);
  printf ("empiric variance:  %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (empiric_mean) < 35);
  TEST_EXPECT (fabsl (empiric_var - sigma * sigma) < 100000);

  fmpz_mod_poly_clear (poly, ctx);
  fmpz_mod_ctx_clear (ctx);
  fmpz_clear (coeff);
  fmpz_clear (halfq);
  fmpz_clear (q);
  fmpz_clear (mod);
}

static void
single_prec (void)
{
  const long double sigma = (long double)1.55 * (1 << SP_LOG2O);
  uint8_t seed[32] = { 0 };
  uint64_t dom = 0;
  uint64_t *seed_ptr = (uint64_t *)&(seed[0]);
  long sum = 0, sum_sqr = 0;
  long double empiric_mean, empiric_var;
  INT_T (sample, 1);
  INT_T (sample_sqr, 2);
  INTVEC_T (samplevec, SAMPLEVEC_DIM, 1);
  unsigned int i;
  int64_t tmp;

  bytes_urandom (seed, sizeof (seed));

  printf ("single precision test\n\n");

  for (*seed_ptr = 0; *seed_ptr < SP_NSAMPLES; (*seed_ptr)++)
    {
      dom = 0;
      int_grandom (sample, SP_LOG2O, seed, dom);
      dom = 1;
      intvec_grandom (samplevec, SP_LOG2O, seed, dom);

      tmp = int_get_i64 (sample);
      sum += tmp;
      sum_sqr += tmp * tmp;

      _VEC_FOREACH_ELEM (samplevec, i)
      {
        tmp = intvec_get_elem_i64 (samplevec, i);
        sum += tmp;
        sum_sqr += tmp * tmp;
      }
    }

  empiric_mean = (long double)sum / (SP_NSAMPLES * (1 + SAMPLEVEC_DIM));
  empiric_var = (long double)sum_sqr / (SP_NSAMPLES * (1 + SAMPLEVEC_DIM));

  printf ("expected mean:     0\n");
  printf ("empiric mean:      %Lf\n", empiric_mean);

  printf ("expected variance: %Lf\n", sigma * sigma);
  printf ("empiric variance:  %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (empiric_mean) < 1);
  TEST_EXPECT (fabsl (empiric_var - sigma * sigma) < 10000);
}

static void
multiple_prec (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  uint64_t *seed_ptr = (uint64_t *)&(seed[0]);
  INT_T (sample, 4);
  INT_T (sum, 4);
  INT_T (sample_sqr, 8);
  INT_T (sum_sqr, 8);

  bytes_urandom (seed, sizeof (seed));

  printf ("multiple precision test\n\n");

  int_set_i64 (sum, 0);
  int_set_i64 (sum_sqr, 0);

  for (*seed_ptr = 0; *seed_ptr < (1 << MP_LOG2NSAMPLES); (*seed_ptr)++)
    {
      int_grandom (sample, MP_LOG2O, seed, dom);

      int_mul (sample_sqr, sample, sample);
      int_add (sum, sum, sample);
      int_add (sum_sqr, sum_sqr, sample_sqr);
    }

  /* divide by number of samples */
  int_rshift (sum, sum, MP_LOG2NSAMPLES);
  int_rshift (sum_sqr, sum_sqr, MP_LOG2NSAMPLES);

  printf ("expected mean:     0\n");
  printf ("empiric mean:      ");
  int_out_str (stdout, 10, sum);
  printf ("\n");

  printf ("expected variance: "
          "2781904943926521595051292914833726986174811381592014551047968455790"
          "11293959946.24\n");
  printf ("empiric variance:  ");
  int_out_str (stdout, 10, sum_sqr);
  printf ("\n\n");
}
