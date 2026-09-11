#include "abdlop-params1.h"
#include "test.h"
#include "src/flint/abdlop.h"
#include <math.h>

#define NSAMPLES 1000000
#define NROWS 3
#define NCOLS 3

#define FLINT_NSAMPLES 500
#define FLINT_NELEMS 16

static void flint_test (void);

int
main (void)
{
  uint8_t seed[32] = { 0 };
  unsigned int i, j, o, k, l, n;
  int64_t tmp;
  uint32_t dom = 0;
  long double mean, empiric_mean;
  long double var, empiric_var;

  lazer_init();

  bytes_urandom (seed, sizeof(seed));

  for (k = 1; k <= 3; k++)
    {
      for (i = 1; i <= NROWS; i++)
        {
          for (j = 1; j <= NCOLS; j++)
            {
              INT_T (z, 2);
              INTVEC_T (v, j, 2);
              INTMAT_T (m, i, j, 2);
              mean = 0;
              var = 0;

              for (o = 0; o < NSAMPLES; o++)
                {
                  int_brandom (z, k, seed, dom);
                  dom++;
                  tmp = int_get_i64 (z);
                  mean += tmp;
                  var += (tmp * tmp);

                  intvec_brandom (v, k, seed, dom);
                  dom++;
                  _VEC_FOREACH_ELEM (v, l)
                  {
                    tmp = intvec_get_elem_i64 (v, l);
                    mean += tmp;
                    var += (tmp * tmp);
                  }

                  intmat_brandom (m, k, seed, dom);
                  dom++;
                  _MAT_FOREACH_ELEM (m, l, n)
                  {
                    tmp = intmat_get_elem_i64 (m, l, n);
                    mean += tmp;
                    var += (tmp * tmp);
                  }
                }

              empiric_mean = mean / (NSAMPLES * (i * j + j + 1));
              empiric_var = var / (NSAMPLES * (i * j + j + 1));

              printf ("bin_%d\n", k);
              printf ("expected mean:     0\n");
              printf ("empiric mean:      %Lf\n", empiric_mean);
              printf ("expected variance: %.2f\n", (float)k / 2);
              printf ("empiric variance:  %Lf\n", empiric_var);
              printf ("\n");

              TEST_EXPECT (fabsl (empiric_mean) < 0.01);
              TEST_EXPECT (fabsl (empiric_var - (float)k / 2) < 0.01);
            }
        }
    }

  flint_test ();

  TEST_PASS ();
}

/* flint test: fmpz_mod_poly_urand_autostable bounded sampling */
static void
flint_test (void)
{
  polyring_srcptr Rq = params1->ring;
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  fmpz_t mod, halfq, c;
  fmpz_mod_ctx_t ctx;
  fmpz_mod_poly_t poly;
  unsigned int it, k;
  const int64_t bnd = 1;
  const unsigned int log2 = 8;
  long double sum = 0, sumsq = 0;

  bytes_urandom (seed, sizeof (seed));

  fmpz_init_int (mod, Rq->q);
  fmpz_init (halfq);
  fmpz_fdiv_q_2exp (halfq, mod, 1);
  fmpz_init (c);

  fmpz_mod_ctx_init (ctx, mod);
  fmpz_mod_poly_init (poly, ctx);

  for (it = 0; it < FLINT_NSAMPLES; it++)
    {
      fmpz_mod_poly_urand_autostable (poly, FLINT_NELEMS, ctx, bnd, log2,
                                      seed, dom);
      dom++;
      for (k = 0; k < FLINT_NELEMS; k++)
        {
          fmpz_mod_poly_get_coeff_fmpz (c, poly, k, ctx);
          if (fmpz_cmp (c, halfq) > 0)
            fmpz_sub (c, c, mod); /* center-reduce to signed */
          {
            int64_t v = fmpz_get_si (c);
            sum += v;
            sumsq += (long double)v * v;
            TEST_EXPECT (v >= -(2 * bnd) && v <= (2 * bnd));
          }
        }
    }

  fmpz_mod_poly_clear (poly, ctx);
  fmpz_mod_ctx_clear (ctx);
  fmpz_clear (c);
  fmpz_clear (halfq);
  fmpz_clear (mod);

  {
    long double n = (long double)FLINT_NSAMPLES * FLINT_NELEMS;
    long double em = sum / n;
    long double ev = sumsq / n;
    printf ("flint bounded test (bnd=%ld)\n", (long)bnd);
    printf ("empiric mean:     %Lf\n", em);
    printf ("empiric variance: %Lf\n", ev);
    printf ("\n");
    TEST_EXPECT (fabsl (em) < 1);
  }
}
