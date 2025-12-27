//! Berlekamp's algorithm for polynomial factorization over finite fields.
//!
//! Factors a squarefree polynomial over GF(p) by computing the Berlekamp
//! subalgebra and splitting with random linear combinations.

use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use std::ops::Neg;
use tertius_poly::dense::DensePoly;
use tertius_rings::finite_field::FiniteField;
use tertius_rings::traits::{Field, Ring};

// Helper functions to avoid trait ambiguity
fn ff_zero<const P: u64>() -> FiniteField<P> {
    <FiniteField<P> as Ring>::zero()
}

fn ff_one<const P: u64>() -> FiniteField<P> {
    <FiniteField<P> as Ring>::one()
}

fn ff_is_zero<const P: u64>(x: &FiniteField<P>) -> bool {
    <FiniteField<P> as Zero>::is_zero(x)
}

/// Result of Berlekamp factorization.
#[derive(Clone, Debug)]
pub struct BerlekampResult<const P: u64> {
    /// Irreducible factors of the input polynomial.
    pub factors: Vec<DensePoly<FiniteField<P>>>,
    /// Number of GCD computations performed.
    pub gcd_count: usize,
}

/// Factors a squarefree polynomial over GF(p) using Berlekamp's algorithm.
///
/// # Arguments
///
/// * `f` - A squarefree polynomial over GF(p)
///
/// # Returns
///
/// The complete factorization into irreducible factors.
///
/// # Panics
///
/// Panics if the polynomial is zero.
pub fn berlekamp_factor<const P: u64>(f: &DensePoly<FiniteField<P>>) -> BerlekampResult<P> {
    assert!(!f.is_zero(), "Cannot factor zero polynomial");

    let n = f.degree();
    if n <= 0 {
        return BerlekampResult {
            factors: vec![f.clone()],
            gcd_count: 0,
        };
    }

    // Compute Berlekamp matrix Q
    let q_matrix = compute_berlekamp_matrix(f);

    // Find null space of (Q - I)
    let null_space = find_null_space(&q_matrix, n);

    // Number of irreducible factors equals dimension of null space
    let num_factors = null_space.len();

    if num_factors == 1 {
        // Polynomial is irreducible
        return BerlekampResult {
            factors: vec![f.clone()],
            gcd_count: 0,
        };
    }

    // Split using null space vectors
    let mut factors = vec![f.clone()];
    let mut gcd_count = 0;
    let mut rng = ChaCha8Rng::seed_from_u64(42);

    while factors.len() < num_factors {
        let mut new_factors = Vec::new();

        for factor in factors {
            if factor.degree() <= 1 {
                new_factors.push(factor);
                continue;
            }

            // Try to split this factor
            if let Some((g1, g2)) = try_split(&factor, &null_space, &mut rng, &mut gcd_count) {
                new_factors.push(g1);
                new_factors.push(g2);
            } else {
                new_factors.push(factor);
            }
        }

        factors = new_factors;
    }

    // Make factors monic
    let factors: Vec<_> = factors.into_iter().map(|f| make_monic(&f)).collect();

    BerlekampResult { factors, gcd_count }
}

/// Computes the Berlekamp matrix Q.
///
/// Q[i,j] = coefficient of x^j in x^(i*p) mod f(x)
fn compute_berlekamp_matrix<const P: u64>(
    f: &DensePoly<FiniteField<P>>,
) -> Vec<Vec<FiniteField<P>>> {
    let n = f.degree();

    // Compute x^p mod f(x) first
    let x = DensePoly::new(vec![ff_zero(), ff_one()]);
    let x_p = pow_mod(&x, P, f);

    // Build matrix: row i contains coefficients of x^(ip) mod f
    let mut matrix = Vec::with_capacity(n);

    // Row 0: x^0 = 1
    let mut row0 = vec![ff_zero(); n];
    row0[0] = ff_one();
    matrix.push(row0);

    // Row i: x^(ip) = (x^p)^i mod f
    let mut current = x_p.clone();
    for _ in 1..n {
        let mut row = vec![ff_zero(); n];
        for i in 0..n.min(current.coeffs().len()) {
            row[i] = current.coeffs()[i].clone();
        }
        matrix.push(row);
        current = poly_mod(&poly_mul(&current, &x_p), f);
    }

    matrix
}

/// Finds the null space of (Q - I).
fn find_null_space<const P: u64>(
    q_matrix: &[Vec<FiniteField<P>>],
    n: usize,
) -> Vec<Vec<FiniteField<P>>> {
    // Compute Q - I
    let mut matrix: Vec<Vec<FiniteField<P>>> = q_matrix.to_vec();
    for i in 0..n {
        matrix[i][i] = matrix[i][i].clone() + ff_one().neg();
    }

    // Gaussian elimination to find null space
    // We want vectors v such that (Q - I)^T v = 0, i.e., v^T (Q - I) = 0
    // Transpose and find null space
    let mut transposed = vec![vec![ff_zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            transposed[i][j] = matrix[j][i].clone();
        }
    }

    // Row reduce the transposed matrix augmented with identity
    let mut augmented = vec![vec![ff_zero(); 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            augmented[i][j] = transposed[i][j].clone();
        }
        augmented[i][n + i] = ff_one();
    }

    // Gaussian elimination
    let mut pivot_col = 0;
    let mut pivot_row = 0;

    while pivot_row < n && pivot_col < n {
        // Find pivot
        let mut max_row = None;
        for row in pivot_row..n {
            if !ff_is_zero(&augmented[row][pivot_col]) {
                max_row = Some(row);
                break;
            }
        }

        if let Some(max_row) = max_row {
            // Swap rows
            augmented.swap(pivot_row, max_row);

            // Scale pivot row
            let pivot_val = augmented[pivot_row][pivot_col].clone();
            let pivot_inv = pivot_val.inv().expect("pivot must be invertible");
            for j in 0..2 * n {
                augmented[pivot_row][j] = augmented[pivot_row][j].clone() * pivot_inv.clone();
            }

            // Eliminate column
            for row in 0..n {
                if row != pivot_row && !ff_is_zero(&augmented[row][pivot_col]) {
                    let factor = augmented[row][pivot_col].clone();
                    for j in 0..2 * n {
                        let sub = factor.clone() * augmented[pivot_row][j].clone();
                        augmented[row][j] = augmented[row][j].clone() + sub.neg();
                    }
                }
            }

            pivot_row += 1;
        }
        pivot_col += 1;
    }

    // Extract null space vectors (rows where left side is all zeros)
    let mut null_space = Vec::new();
    for i in 0..n {
        let is_zero_row = (0..n).all(|j| ff_is_zero(&augmented[i][j]));
        if is_zero_row {
            let vec: Vec<_> = (n..2 * n).map(|j| augmented[i][j].clone()).collect();
            null_space.push(vec);
        }
    }

    // Always include the constant vector (1, 0, 0, ..., 0)
    if null_space.is_empty() {
        let mut v = vec![ff_zero(); n];
        v[0] = ff_one();
        null_space.push(v);
    }

    null_space
}

/// Tries to split a polynomial using null space vectors.
fn try_split<const P: u64, R: Rng>(
    f: &DensePoly<FiniteField<P>>,
    null_space: &[Vec<FiniteField<P>>],
    rng: &mut R,
    gcd_count: &mut usize,
) -> Option<(DensePoly<FiniteField<P>>, DensePoly<FiniteField<P>>)> {
    // Try random linear combinations of null space vectors
    for _ in 0..50 {
        // Build random linear combination
        let mut coeffs = vec![ff_zero(); f.degree()];
        for vec in null_space {
            let rand_coeff = FiniteField::new(rng.gen_range(0..P));
            for (i, v) in vec.iter().enumerate() {
                if i < coeffs.len() {
                    coeffs[i] = coeffs[i].clone() + rand_coeff.clone() * v.clone();
                }
            }
        }

        let h = DensePoly::new(coeffs);
        if h.is_zero() || h.degree() <= 0 {
            continue;
        }

        // Try h^((p-1)/2) - 1 for splitting (works for odd p)
        if P > 2 {
            let exp = (P - 1) / 2;
            let h_pow = pow_mod(&h, exp, f);
            let h_pow_minus_1 = poly_sub(&h_pow, &DensePoly::new(vec![ff_one()]));

            *gcd_count += 1;
            let g = poly_gcd(f, &h_pow_minus_1);

            if g.degree() > 0 && g.degree() < f.degree() {
                let other = poly_div(f, &g).0;
                return Some((g, other));
            }

            // Also try h^((p-1)/2) + 1
            let h_pow_plus_1 = poly_add(&h_pow, &DensePoly::new(vec![ff_one()]));
            *gcd_count += 1;
            let g = poly_gcd(f, &h_pow_plus_1);

            if g.degree() > 0 && g.degree() < f.degree() {
                let other = poly_div(f, &g).0;
                return Some((g, other));
            }
        } else {
            // For p = 2, just use GCD with h
            *gcd_count += 1;
            let g = poly_gcd(f, &h);
            if g.degree() > 0 && g.degree() < f.degree() {
                let other = poly_div(f, &g).0;
                return Some((g, other));
            }
        }
    }

    None
}

/// Computes a^exp mod f using binary exponentiation.
fn pow_mod<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    mut exp: u64,
    f: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    if exp == 0 {
        return DensePoly::new(vec![ff_one()]);
    }

    let mut result = DensePoly::new(vec![ff_one()]);
    let mut base = poly_mod(a, f);

    while exp > 0 {
        if exp & 1 == 1 {
            result = poly_mod(&poly_mul(&result, &base), f);
        }
        base = poly_mod(&poly_mul(&base, &base), f);
        exp >>= 1;
    }

    result
}

/// Polynomial multiplication.
fn poly_mul<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    b: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![ff_zero(); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        for (j, bj) in b_coeffs.iter().enumerate() {
            result[i + j] = result[i + j].clone() + ai.clone() * bj.clone();
        }
    }

    DensePoly::new(result)
}

/// Polynomial addition.
fn poly_add<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    b: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![ff_zero(); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i] = c.clone();
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        result[i] = result[i].clone() + c.clone();
    }

    DensePoly::new(result)
}

/// Polynomial subtraction.
fn poly_sub<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    b: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![ff_zero(); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i] = c.clone();
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        result[i] = result[i].clone() + c.clone().neg();
    }

    DensePoly::new(result)
}

/// Polynomial modular reduction.
fn poly_mod<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    f: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    poly_div(a, f).1
}

/// Polynomial division with remainder.
fn poly_div<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    b: &DensePoly<FiniteField<P>>,
) -> (DensePoly<FiniteField<P>>, DensePoly<FiniteField<P>>) {
    if b.is_zero() {
        panic!("Division by zero polynomial");
    }

    if a.degree() < b.degree() {
        return (DensePoly::zero(), a.clone());
    }

    let mut remainder = a.coeffs().to_vec();
    let divisor = b.coeffs();
    let divisor_lead_inv = divisor.last().unwrap().inv().expect("divisor lead must be invertible");
    let deg_diff = a.degree() - b.degree();
    let mut quotient = vec![ff_zero(); deg_diff + 1];

    for i in (0..=deg_diff).rev() {
        if remainder.len() > i + divisor.len() - 1 {
            let coeff = remainder[i + divisor.len() - 1].clone() * divisor_lead_inv.clone();
            quotient[i] = coeff.clone();

            for (j, d) in divisor.iter().enumerate() {
                let idx = i + j;
                if idx < remainder.len() {
                    remainder[idx] = remainder[idx].clone() + (coeff.clone() * d.clone()).neg();
                }
            }
        }
    }

    // Trim trailing zeros
    while remainder.len() > 1 && remainder.last().map_or(false, |c| ff_is_zero(c)) {
        remainder.pop();
    }

    (DensePoly::new(quotient), DensePoly::new(remainder))
}

/// Polynomial GCD using Euclidean algorithm.
fn poly_gcd<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    b: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    let mut a = a.clone();
    let mut b = b.clone();

    while !b.is_zero() {
        let r = poly_div(&a, &b).1;
        a = b;
        b = r;
    }

    make_monic(&a)
}

/// Makes a polynomial monic (leading coefficient = 1).
fn make_monic<const P: u64>(f: &DensePoly<FiniteField<P>>) -> DensePoly<FiniteField<P>> {
    if f.is_zero() {
        return f.clone();
    }

    let coeffs = f.coeffs();
    let lead_inv = coeffs.last().unwrap().inv().expect("lead must be invertible");
    let new_coeffs: Vec<_> = coeffs.iter().map(|c| c.clone() * lead_inv.clone()).collect();

    DensePoly::new(new_coeffs)
}

/// Parallel factorization of multiple polynomials.
pub fn berlekamp_factor_batch<const P: u64>(
    polys: &[DensePoly<FiniteField<P>>],
) -> Vec<BerlekampResult<P>> {
    polys.par_iter().map(berlekamp_factor).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    type GF5 = FiniteField<5>;
    type GF7 = FiniteField<7>;

    fn gf5(n: i64) -> GF5 {
        GF5::new(n.rem_euclid(5) as u64)
    }

    fn gf7(n: i64) -> GF7 {
        GF7::new(n.rem_euclid(7) as u64)
    }

    #[test]
    fn test_poly_mul() {
        // (x + 1)(x + 2) = x^2 + 3x + 2 in GF(5)
        let a = DensePoly::new(vec![gf5(1), gf5(1)]); // x + 1
        let b = DensePoly::new(vec![gf5(2), gf5(1)]); // x + 2
        let c = poly_mul(&a, &b);

        let expected = DensePoly::new(vec![gf5(2), gf5(3), gf5(1)]);
        assert_eq!(c.coeffs(), expected.coeffs());
    }

    #[test]
    fn test_poly_gcd() {
        // GCD of (x+1)^2 and (x+1)(x+2) should be (x+1)
        let a = poly_mul(
            &DensePoly::new(vec![gf5(1), gf5(1)]),
            &DensePoly::new(vec![gf5(1), gf5(1)]),
        );
        let b = poly_mul(
            &DensePoly::new(vec![gf5(1), gf5(1)]),
            &DensePoly::new(vec![gf5(2), gf5(1)]),
        );

        let g = poly_gcd(&a, &b);
        assert_eq!(g.degree(), 1);
    }

    #[test]
    fn test_factor_linear() {
        // x + 1 is irreducible
        let f = DensePoly::new(vec![gf5(1), gf5(1)]);
        let result = berlekamp_factor(&f);
        assert_eq!(result.factors.len(), 1);
    }

    #[test]
    fn test_factor_reducible() {
        // (x + 1)(x + 2) = x^2 + 3x + 2 over GF(5)
        let f = DensePoly::new(vec![gf5(2), gf5(3), gf5(1)]);
        let result = berlekamp_factor(&f);

        // Should factor into two linear factors
        assert_eq!(result.factors.len(), 2);

        // Verify product equals original
        let product = poly_mul(&result.factors[0], &result.factors[1]);
        let f_monic = make_monic(&f);
        assert_eq!(product.coeffs(), f_monic.coeffs());
    }

    #[test]
    fn test_factor_irreducible_quadratic() {
        // x^2 + 2 is irreducible over GF(5) (2 is not a quadratic residue mod 5)
        let f = DensePoly::new(vec![gf5(2), gf5(0), gf5(1)]);
        let result = berlekamp_factor(&f);
        assert_eq!(result.factors.len(), 1);
    }

    #[test]
    fn test_factor_three_factors() {
        // (x)(x + 1)(x + 2) = x^3 + 3x^2 + 2x over GF(5)
        let x = DensePoly::new(vec![gf5(0), gf5(1)]);
        let xp1 = DensePoly::new(vec![gf5(1), gf5(1)]);
        let xp2 = DensePoly::new(vec![gf5(2), gf5(1)]);

        let f = poly_mul(&poly_mul(&x, &xp1), &xp2);
        let result = berlekamp_factor(&f);

        assert_eq!(result.factors.len(), 3);
    }

    #[test]
    fn test_pow_mod() {
        let f = DensePoly::new(vec![gf7(1), gf7(0), gf7(1)]); // x^2 + 1
        let x = DensePoly::new(vec![gf7(0), gf7(1)]); // x

        // x^7 mod (x^2 + 1) over GF(7)
        let result = pow_mod(&x, 7, &f);
        assert!(result.degree() < 2);
    }
}
