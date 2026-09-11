#include "test.h"
#include "src/flint/abdlop.h"
#include <math.h>

#define NSAMPLES 1000000
#define SAMPLEVEC_DIM 4

static inline long double
_expected_mean (int64_t lo, int64_t hi)
{
  return ((long double)lo + hi) / 2;
}

static inline long double
_expected_var (int64_t lo, int64_t hi)
{
  return (long double)(((hi - lo + 1) * (hi - lo + 1)) - 1) / 12;
}

static void test_urandom_bnd (int64_t lo, int64_t hi);
static void test_urandom (int64_t mod, unsigned int log2mod);
static void test_urandom_flint (int64_t mod, unsigned int log2mod);

int
main (void)
{
  lazer_init();

  test_urandom (3, 2);
  test_urandom (15, 4);
  test_urandom (29, 5);

  test_urandom_flint (3, 2);
  test_urandom_flint (15, 4);
  test_urandom_flint (29, 5);

  test_urandom_bnd (-2, 2);
  test_urandom_bnd (-5, -4);
  test_urandom_bnd (15, 16);
  test_urandom_bnd (-5, 4);
  test_urandom_bnd (-4, 5);
  test_urandom_bnd (-4, 0);
  test_urandom_bnd (0, 4);

  TEST_PASS ();
}

/* flint test: fmpz_mod_poly_urand must match lazer intvec_urandom */
static void
test_urandom_flint (int64_t mod, unsigned int log2mod)
{
  unsigned int i;
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  INT_T (mod_i, 1);
  INTVEC_T (v, SAMPLEVEC_DIM, 1);
  fmpz_t mod_f;
  fmpz_mod_ctx_t ctx;
  fmpz_mod_poly_t poly;
  unsigned int k;

  bytes_urandom (seed, sizeof (seed));
  int_set_i64 (mod_i, mod);
  fmpz_init_set_ui (mod_f, mod);
  fmpz_mod_ctx_init (ctx, mod_f);
  fmpz_mod_poly_init (poly, ctx);

  for (i = 0; i < 100; i++)
    {
      intvec_urandom (v, mod_i, log2mod, seed, dom);
      fmpz_mod_poly_urand (poly, SAMPLEVEC_DIM, ctx, mod_f, log2mod, seed,
                           dom);

      for (k = 0; k < SAMPLEVEC_DIM; k++)
        {
          fmpz_t c;
          fmpz_init (c);
          fmpz_mod_poly_get_coeff_fmpz (c, poly, k, ctx);
          TEST_EXPECT (fmpz_get_ui (c) == (ulong)intvec_get_elem_i64 (v, k));
          fmpz_clear (c);
        }
      dom++;
    }

  fmpz_mod_poly_clear (poly, ctx);
  fmpz_mod_ctx_clear (ctx);
  fmpz_clear (mod_f);
}

static void
test_urandom_bnd (int64_t lo, int64_t hi)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  uint64_t *seed_ptr = (uint64_t *)&(seed[0]);
  int64_t z;
  long double sum, sum_sqr, empiric_mean, empiric_var, expected_mean,
      expected_var;
  INT_T (low, 3);
  INT_T (high, 3);
  INT_T (sample, 3);
  INTVEC_T (samplevec, SAMPLEVEC_DIM, 3);
  unsigned int i;

  bytes_urandom (seed, sizeof (seed));

  sum = 0;
  sum_sqr = 0;
  int_set_i64 (low, lo);
  int_set_i64 (high, hi);

  expected_mean = _expected_mean (lo, hi);
  expected_var = _expected_var (lo, hi);

  for (*seed_ptr = 0; *seed_ptr < NSAMPLES; (*seed_ptr)++)
    {
      dom = 0;
      int_urandom_bnd (sample, low, high, seed, dom);
      z = int_get_i64 (sample);
      TEST_EXPECT (sample->limbs[1] == 0);
      TEST_EXPECT (sample->limbs[2] == 0);
      TEST_EXPECT (z <= hi);
      TEST_EXPECT (z >= lo);

      sum += z;
      sum_sqr += ((long double)z - expected_mean)
                 * ((long double)z - expected_mean);

      dom = 1;
      intvec_urandom_bnd (samplevec, low, high, seed, dom);
      _VEC_FOREACH_ELEM (samplevec, i)
      {
        z = intvec_get_elem_i64 (samplevec, i);
        TEST_EXPECT (z <= hi);
        TEST_EXPECT (z >= lo);

        sum += z;
        sum_sqr += ((long double)z - expected_mean)
                   * ((long double)z - expected_mean);
      }
    }

  empiric_mean = sum / (NSAMPLES * (1 + SAMPLEVEC_DIM));
  empiric_var = sum_sqr / (NSAMPLES * (1 + SAMPLEVEC_DIM));

  printf ("uniform in [%ld,%ld]\n", lo, hi);
  printf ("expected mean:     %Lf\n", expected_mean);
  printf ("empiric mean:      %Lf\n", empiric_mean);
  printf ("expected variance: %Lf\n", expected_var);
  printf ("empiric variance:  %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (expected_mean - empiric_mean) < 0.01);
  TEST_EXPECT (fabsl (expected_var - empiric_var) < 0.01);
}

static void
test_urandom (int64_t mod, unsigned int log2mod)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  uint64_t *seed_ptr = (uint64_t *)&(seed[0]);
  int64_t z;
  long double sum, sum_sqr, empiric_mean, empiric_var, expected_mean,
      expected_var;
  INT_T (tmp, 1);
  INT_T (sample, 3);
  INTVEC_T (samplevec, SAMPLEVEC_DIM, 3);
  unsigned int i;

  bytes_urandom (seed, sizeof (seed));
  int_set_i64 (tmp, mod);

  sum = 0;
  sum_sqr = 0;

  expected_mean = _expected_mean (0, mod - 1);
  expected_var = _expected_var (0, mod - 1);

  for (*seed_ptr = 0; *seed_ptr < NSAMPLES; (*seed_ptr)++)
    {
      dom = 0;
      int_urandom (sample, tmp, log2mod, seed, dom);
      z = int_get_i64 (sample);
      TEST_EXPECT (sample->limbs[1] == 0);
      TEST_EXPECT (sample->limbs[2] == 0);
      TEST_EXPECT (z <= (mod - 1));

      sum += z;
      sum_sqr += ((long double)z - expected_mean)
                 * ((long double)z - expected_mean);

      dom = 1;
      intvec_urandom (samplevec, tmp, log2mod, seed, dom);
      _VEC_FOREACH_ELEM (samplevec, i)
      {
        z = intvec_get_elem_i64 (samplevec, i);
        TEST_EXPECT (z <= (mod - 1));
        sum += z;
        sum_sqr += ((long double)z - expected_mean)
                   * ((long double)z - expected_mean);
      }
    }

  empiric_mean = sum / (NSAMPLES * (1 + SAMPLEVEC_DIM));
  empiric_var = sum_sqr / (NSAMPLES * (1 + SAMPLEVEC_DIM));

  printf ("uniform in [0,%ld]\n", mod - 1);
  printf ("expected mean:     %Lf\n", expected_mean);
  printf ("empiric mean:      %Lf\n", empiric_mean);
  printf ("expected variance: %Lf\n", expected_var);
  printf ("empiric variance:  %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (expected_mean - empiric_mean) < 0.1);
  TEST_EXPECT (fabsl (expected_var - empiric_var) < 0.1);
}
