//! Cantor-Zassenhaus algorithm for polynomial factorization over finite fields.
//!
//! This is a probabilistic algorithm for equal-degree factorization that
//! works well for large primes. For complete factorization, use in conjunction
//! with distinct-degree factorization.

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

/// Result of Cantor-Zassenhaus factorization.
#[derive(Clone, Debug)]
pub struct CantorZassenhausResult<const P: u64> {
    /// Irreducible factors (all of the same degree if equal-degree factorization).
    pub factors: Vec<DensePoly<FiniteField<P>>>,
    /// Number of random attempts made.
    pub attempts: usize,
}

/// Performs complete factorization using distinct-degree then equal-degree factorization.
pub fn cantor_zassenhaus_factor<const P: u64>(
    f: &DensePoly<FiniteField<P>>,
) -> CantorZassenhausResult<P> {
    if f.is_zero() {
        return CantorZassenhausResult {
            factors: vec![],
            attempts: 0,
        };
    }

    let n = f.degree();
    if n == 0 {
        return CantorZassenhausResult {
            factors: vec![f.clone()],
            attempts: 0,
        };
    }

    if n == 1 {
        return CantorZassenhausResult {
            factors: vec![make_monic(f)],
            attempts: 0,
        };
    }

    // Step 1: Distinct-degree factorization
    let ddf_result = distinct_degree_factorization(f);

    // Step 2: Equal-degree factorization for each component
    let mut all_factors = Vec::new();
    let mut total_attempts = 0;

    for (degree, poly) in ddf_result {
        if poly.degree() == degree {
            all_factors.push(make_monic(&poly));
        } else if poly.degree() > 0 {
            let edf_result = equal_degree_factorization(&poly, degree);
            total_attempts += edf_result.attempts;
            all_factors.extend(edf_result.factors);
        }
    }

    CantorZassenhausResult {
        factors: all_factors,
        attempts: total_attempts,
    }
}

/// Distinct-degree factorization.
fn distinct_degree_factorization<const P: u64>(
    f: &DensePoly<FiniteField<P>>,
) -> Vec<(usize, DensePoly<FiniteField<P>>)> {
    let n = f.degree();
    let mut result = Vec::new();
    let mut h = f.clone();

    let x = DensePoly::new(vec![ff_zero(), ff_one()]);
    let mut x_pow = pow_mod(&x, P, f);

    for d in 1..=n / 2 {
        if h.degree() < 2 * d {
            break;
        }

        if d > 1 {
            x_pow = pow_mod(&x_pow, P, &h);
        }

        let x_pow_minus_x = poly_sub(&x_pow, &x);
        let g = poly_gcd(&h, &x_pow_minus_x);

        if g.degree() > 0 {
            result.push((d, g.clone()));
            h = poly_div(&h, &g).0;
            x_pow = poly_mod(&x_pow, &h);
        }
    }

    if h.degree() > 0 {
        result.push((h.degree(), h));
    }

    result
}

/// Equal-degree factorization using Cantor-Zassenhaus.
fn equal_degree_factorization<const P: u64>(
    f: &DensePoly<FiniteField<P>>,
    d: usize,
) -> CantorZassenhausResult<P> {
    let n = f.degree();

    if n == d {
        return CantorZassenhausResult {
            factors: vec![make_monic(f)],
            attempts: 0,
        };
    }

    let num_factors = n / d;
    let mut factors = vec![f.clone()];
    let mut attempts = 0;
    let mut rng = ChaCha8Rng::seed_from_u64(123);

    while factors.len() < num_factors {
        let mut new_factors = Vec::new();

        for factor in &factors {
            if factor.degree() == d {
                new_factors.push(factor.clone());
                continue;
            }

            if let Some((g1, g2)) = try_split_equal_degree(factor, d, &mut rng, &mut attempts) {
                new_factors.push(g1);
                new_factors.push(g2);
            } else {
                new_factors.push(factor.clone());
            }
        }

        factors = new_factors;

        if attempts > 1000 {
            break;
        }
    }

    let factors = factors.into_iter().map(|f| make_monic(&f)).collect();

    CantorZassenhausResult { factors, attempts }
}

fn try_split_equal_degree<const P: u64, R: Rng>(
    f: &DensePoly<FiniteField<P>>,
    _d: usize,
    rng: &mut R,
    attempts: &mut usize,
) -> Option<(DensePoly<FiniteField<P>>, DensePoly<FiniteField<P>>)> {
    let n = f.degree();

    for _ in 0..20 {
        *attempts += 1;

        let mut coeffs = Vec::with_capacity(n);
        for _ in 0..n {
            coeffs.push(FiniteField::new(rng.gen_range(0..P)));
        }
        let a = DensePoly::new(coeffs);

        if a.is_zero() || a.degree() == 0 {
            continue;
        }

        let g = poly_gcd(&a, f);
        if g.degree() > 0 && g.degree() < f.degree() {
            let other = poly_div(f, &g).0;
            return Some((g, other));
        }

        if P > 2 {
            let exp = (P - 1) / 2;
            let a_pow = pow_mod(&a, exp, f);

            let a_pow_minus_1 = poly_sub(&a_pow, &DensePoly::new(vec![ff_one()]));
            let g = poly_gcd(f, &a_pow_minus_1);

            if g.degree() > 0 && g.degree() < f.degree() {
                let other = poly_div(f, &g).0;
                return Some((g, other));
            }

            let a_pow_plus_1 = poly_add(&a_pow, &DensePoly::new(vec![ff_one()]));
            let g = poly_gcd(f, &a_pow_plus_1);

            if g.degree() > 0 && g.degree() < f.degree() {
                let other = poly_div(f, &g).0;
                return Some((g, other));
            }
        } else {
            let mut trace = a.clone();
            let mut current = a.clone();
            for _ in 1..n {
                current = poly_mod(&poly_mul(&current, &current), f);
                trace = poly_add(&trace, &current);
            }

            let g = poly_gcd(f, &trace);
            if g.degree() > 0 && g.degree() < f.degree() {
                let other = poly_div(f, &g).0;
                return Some((g, other));
            }
        }
    }

    None
}

pub fn cantor_zassenhaus_factor_batch<const P: u64>(
    polys: &[DensePoly<FiniteField<P>>],
) -> Vec<CantorZassenhausResult<P>> {
    polys.par_iter().map(cantor_zassenhaus_factor).collect()
}

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

fn poly_mod<const P: u64>(
    a: &DensePoly<FiniteField<P>>,
    f: &DensePoly<FiniteField<P>>,
) -> DensePoly<FiniteField<P>> {
    poly_div(a, f).1
}

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
        let idx = i + divisor.len() - 1;
        if idx < remainder.len() {
            let coeff = remainder[idx].clone() * divisor_lead_inv.clone();
            quotient[i] = coeff.clone();

            for (j, d) in divisor.iter().enumerate() {
                let ridx = i + j;
                if ridx < remainder.len() {
                    remainder[ridx] = remainder[ridx].clone() + (coeff.clone() * d.clone()).neg();
                }
            }
        }
    }

    while remainder.len() > 1 && remainder.last().map_or(false, |c| ff_is_zero(c)) {
        remainder.pop();
    }

    (DensePoly::new(quotient), DensePoly::new(remainder))
}

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

fn make_monic<const P: u64>(f: &DensePoly<FiniteField<P>>) -> DensePoly<FiniteField<P>> {
    if f.is_zero() {
        return f.clone();
    }

    let coeffs = f.coeffs();
    let lead_inv = coeffs.last().unwrap().inv().expect("lead must be invertible");
    let new_coeffs: Vec<_> = coeffs.iter().map(|c| c.clone() * lead_inv.clone()).collect();

    DensePoly::new(new_coeffs)
}

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

#[cfg(test)]
mod tests {
    use super::*;

    type GF7 = FiniteField<7>;

    fn gf7(n: i64) -> GF7 {
        GF7::new(n.rem_euclid(7) as u64)
    }

    #[test]
    fn test_factor_linear() {
        let f = DensePoly::new(vec![gf7(1), gf7(1)]);
        let result = cantor_zassenhaus_factor(&f);
        assert_eq!(result.factors.len(), 1);
        assert_eq!(result.factors[0].degree(), 1);
    }

    #[test]
    fn test_factor_two_linears() {
        let f = DensePoly::new(vec![gf7(2), gf7(3), gf7(1)]);
        let result = cantor_zassenhaus_factor(&f);
        assert_eq!(result.factors.len(), 2);

        let product = poly_mul(&result.factors[0], &result.factors[1]);
        assert_eq!(product.coeffs(), f.coeffs());
    }

    #[test]
    fn test_factor_cubic() {
        let x = DensePoly::new(vec![gf7(0), gf7(1)]);
        let xp1 = DensePoly::new(vec![gf7(1), gf7(1)]);
        let xp2 = DensePoly::new(vec![gf7(2), gf7(1)]);

        let f = poly_mul(&poly_mul(&x, &xp1), &xp2);
        let result = cantor_zassenhaus_factor(&f);

        assert_eq!(result.factors.len(), 3);
    }
}
