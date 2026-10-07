#ifndef ARITH_H
#define ARITH_H
/**
 * \file   arith.h
 * \brief  Pluggable arithmetic backend interface.
 *
 * The protocol logic (abdlop / lnp / lin) performs all of its arithmetic
 * through this interface instead of touching backend structs directly.
 * Each backend implements the same op set over the ring
 *
 *     R_q = (Z/qZ)[x] / (x^d + 1),
 *
 * so straight-line protocol code yields identical *values* in every backend.
 *
 * Determinism contract: input generation (sampling, RNG, hashing) lives
 * ABOVE this layer, driven only by ring_new(q,d) and the coefficient
 * accessors.  Backends must not sample internally.
 *
 * Value contract: freshly allocated poly/vec/mat objects hold the zero
 * ring element; *{set,get}_coeff move canonical representatives in the
 * range [0, q); the accumulating ops w <- w + ... are the multiplicative
 * atom any straight-line protocol reduces to, so different backends agree
 * coefficient-wise on the *same* input series.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

  /* Opaque per-backend vtable. */
  typedef struct arith_ops arith_ops_s;
  typedef const arith_ops_s *arith_ops_ptr;

  /* Opaque handles. */
  typedef struct arith_ring arith_ring_s;
  typedef arith_ring_s *arith_ring_ptr;
  typedef struct arith_poly arith_poly_s;
  typedef arith_poly_s *arith_poly_ptr;
  typedef struct arith_vec arith_vec_s;
  typedef arith_vec_s *arith_vec_ptr;
  typedef struct arith_mat arith_mat_s;
  typedef arith_mat_s *arith_mat_ptr;

  struct arith_ops
  {
    const char *name;

    /* ---- ring lifecycle ---- */
    /* spec is optional, backend-specific ring construction data (see
     * backend docs); NULL means "build from q,d alone". */
    arith_ring_ptr (*ring_new) (unsigned long long q, unsigned int d,
                                const void *spec);
    void (*ring_free) (arith_ring_ptr);
    unsigned long long (*ring_q) (const arith_ring_ptr);
    unsigned int (*ring_d) (const arith_ring_ptr);

    /* ---- ring element (a polynomial of degree < d) ---- */
    arith_poly_ptr (*poly_alloc) (const arith_ring_ptr);
    void (*poly_free) (arith_poly_ptr);
    void (*poly_set_coeff) (arith_poly_ptr, unsigned int i, long long c);
    long long (*poly_get_coeff) (const arith_poly_ptr, unsigned int i);

    /* ---- vector of ring elements ---- */
    arith_vec_ptr (*vec_alloc) (const arith_ring_ptr, unsigned int nelems);
    void (*vec_free) (arith_vec_ptr);
    void (*vec_set) (arith_vec_ptr dst, const arith_vec_ptr src);
    void (*vec_set_coeff) (arith_vec_ptr, unsigned int elem,
                           unsigned int k, long long c);
    long long (*vec_get_coeff) (const arith_vec_ptr, unsigned int elem,
                                unsigned int k);

    /* ---- matrix of ring elements ---- */
    arith_mat_ptr (*mat_alloc) (const arith_ring_ptr, unsigned int nrows,
                                unsigned int ncols);
    void (*mat_free) (arith_mat_ptr);
    void (*mat_set_coeff) (arith_mat_ptr, unsigned int i, unsigned int j,
                           unsigned int k, long long c);
    long long (*mat_get_coeff) (const arith_mat_ptr, unsigned int i,
                                unsigned int j, unsigned int k);

    /* ---- linear algebra (accumulating forms) ---- */
    /* w <- w + A * v
     * A: nrows x ncols, v: ncols, w: nrows */
    void (*vec_addmul) (const arith_vec_ptr w, const arith_mat_ptr A,
                        const arith_vec_ptr v);
    /* w <- w + s * v   (s a ring element, coefficient-wise) */
    void (*vec_acc_polyscalar) (const arith_vec_ptr w, const arith_poly_ptr s,
                                const arith_vec_ptr v);
    /* w <- w + a, or w <- w - a when neg != 0 */
    void (*vec_acc) (const arith_vec_ptr w, const arith_vec_ptr a, int neg);
  };

  /* Backend instances. */
  arith_ops_ptr arith_lazer_ops (void);
  arith_ops_ptr arith_flint_ops (void);
  arith_ops_ptr arith_nmod_ops (void); /* FLINT nmod_poly, single limb */
  arith_ops_ptr arith_ntl_ops (void);  /* NTL ZZ_pX */

#ifdef __cplusplus
}
#endif

#endif /* ARITH_H */