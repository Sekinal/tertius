//! Partial fraction decomposition.
//!
//! Given a proper rational function P(x)/Q(x) where deg(P) < deg(Q),
//! decompose into a sum of simpler fractions.
//!
//! If Q(x) = q₁(x)^{e₁} * q₂(x)^{e₂} * ... * qₖ(x)^{eₖ} (factored over the base field),
//! then:
//!
//! P(x)/Q(x) = Σᵢ Σⱼ₌₁^{eᵢ} aᵢⱼ(x) / qᵢ(x)^j
//!
//! where deg(aᵢⱼ) < deg(qᵢ).

use tertius_poly::algorithms::gcd::{poly_div_rem, poly_extended_gcd};
use tertius_poly::dense::DensePoly;
use tertius_rings::traits::Field;

use crate::RationalFunction;

/// A single term in the partial fraction decomposition.
#[derive(Clone, Debug)]
pub struct PartialFractionTerm<K: Field> {
    /// Numerator polynomial (degree < degree of denominator base).
    pub numerator: DensePoly<K>,
    /// The irreducible factor in the denominator.
    pub denominator_base: DensePoly<K>,
    /// The power of the denominator factor.
    pub power: u32,
}

/// Result of partial fraction decomposition.
#[derive(Clone, Debug)]
pub struct PartialFractionDecomposition<K: Field> {
    /// The polynomial part (if deg(num) >= deg(den) before decomposition).
    pub polynomial_part: DensePoly<K>,
    /// The partial fraction terms.
    pub terms: Vec<PartialFractionTerm<K>>,
}

impl<K: Field> PartialFractionDecomposition<K> {
    /// Reconstructs the original rational function from the decomposition.
    pub fn to_rational_function(&self) -> RationalFunction<K> {
        let mut result = RationalFunction::from_poly(self.polynomial_part.clone());

        for term in &self.terms {
            let den = term.denominator_base.pow(term.power);
            let term_rf = RationalFunction::new(term.numerator.clone(), den);
            result = result + term_rf;
        }

        result
    }
}

/// Performs partial fraction decomposition given a factorization of the denominator.
///
/// # Arguments
///
/// * `rf` - The rational function to decompose
/// * `factors` - Factorization of the denominator as (factor, multiplicity) pairs.
///               Factors should be coprime.
///
/// # Returns
///
/// The partial fraction decomposition.
///
/// # Algorithm
///
/// For coprime factors, we use the extended Euclidean algorithm (Bézout identity).
/// For powers, we use successive reduction.
pub fn partial_fraction_decomposition<K: Field>(
    rf: &RationalFunction<K>,
    factors: &[(DensePoly<K>, u32)],
) -> PartialFractionDecomposition<K> {
    // First, extract polynomial part if necessary
    let (poly_part, proper) = rf.decompose_proper();

    if proper.is_zero() {
        return PartialFractionDecomposition {
            polynomial_part: poly_part,
            terms: Vec::new(),
        };
    }

    // Handle single factor case
    if factors.len() == 1 {
        let (factor, power) = &factors[0];
        let terms = decompose_single_power(&proper.numerator().clone(), factor, *power);
        return PartialFractionDecomposition {
            polynomial_part: poly_part,
            terms,
        };
    }

    // Multiple coprime factors: use Bézout identity recursively
    let terms = decompose_coprime_factors(proper.numerator(), proper.denominator(), factors);

    PartialFractionDecomposition {
        polynomial_part: poly_part,
        terms,
    }
}

/// Decomposes A / q^e into A₁/q + A₂/q² + ... + Aₑ/qᵉ
fn decompose_single_power<K: Field>(
    a: &DensePoly<K>,
    q: &DensePoly<K>,
    power: u32,
) -> Vec<PartialFractionTerm<K>> {
    if power == 0 {
        return Vec::new();
    }

    if power == 1 {
        // Already in simplest form - just reduce A mod q
        let (_, remainder) = poly_div_rem(a, q);
        return vec![PartialFractionTerm {
            numerator: remainder,
            denominator_base: q.clone(),
            power: 1,
        }];
    }

    // For higher powers, we decompose iteratively:
    // A / q^e = A₁/q + A₂/q² + ... + Aₑ/qᵉ
    //
    // At each step: A = Q*q + R where deg(R) < deg(q)
    // Then A/qᵉ = Q/qᵉ⁻¹ + R/qᵉ
    // So Aₑ = R, and we continue with Q for powers 1..e-1

    let mut terms = Vec::new();
    let mut current = a.clone();
    let mut current_power = power;

    while current_power > 0 {
        let (quotient, remainder) = poly_div_rem(&current, q);

        terms.push(PartialFractionTerm {
            numerator: remainder,
            denominator_base: q.clone(),
            power: current_power,
        });

        current = quotient;
        current_power -= 1;

        if current.is_zero() {
            break;
        }
    }

    terms
}

/// Decomposes using Bézout identity for coprime factors.
fn decompose_coprime_factors<K: Field>(
    a: &DensePoly<K>,
    _d: &DensePoly<K>,
    factors: &[(DensePoly<K>, u32)],
) -> Vec<PartialFractionTerm<K>> {
    if factors.is_empty() {
        return Vec::new();
    }

    if factors.len() == 1 {
        let (q, e) = &factors[0];
        return decompose_single_power(a, q, *e);
    }

    // Split factors into two groups
    let mid = factors.len() / 2;
    let (left_factors, right_factors) = factors.split_at(mid);

    // Compute products for left and right groups
    let left_product = factor_product(left_factors);
    let right_product = factor_product(right_factors);

    // Use extended GCD: gcd(left, right) = 1 (coprime)
    // So there exist s, t such that s*left + t*right = 1
    // Then A = A*s*left + A*t*right
    // A / (left*right) = A*t / left + A*s / right

    let (gcd, s, t) = poly_extended_gcd(&left_product, &right_product);

    // Verify coprimality (gcd should be constant)
    debug_assert!(gcd.degree() == 0, "Factors must be coprime");

    // Scale by gcd inverse to get proper Bézout coefficients
    let gcd_inv = gcd.coeff(0).inv().expect("gcd should be nonzero");

    // A*t goes with left_product (its denominator is left_product)
    let left_numerator = a.mul(&t).scale(&gcd_inv);
    // A*s goes with right_product
    let right_numerator = a.mul(&s).scale(&gcd_inv);

    // Reduce left_numerator mod left_product and recursively decompose
    let (_, left_reduced) = poly_div_rem(&left_numerator, &left_product);
    let left_terms = decompose_coprime_factors(&left_reduced, &left_product, left_factors);

    // Same for right
    let (_, right_reduced) = poly_div_rem(&right_numerator, &right_product);
    let right_terms = decompose_coprime_factors(&right_reduced, &right_product, right_factors);

    let mut terms = left_terms;
    terms.extend(right_terms);
    terms
}

/// Computes the product of factors with their powers.
fn factor_product<K: Field>(factors: &[(DensePoly<K>, u32)]) -> DensePoly<K> {
    let mut result = DensePoly::one();
    for (f, e) in factors {
        result = result.mul(&f.pow(*e));
    }
    result
}

/// Simple partial fraction for two coprime factors.
///
/// Given A/(B*C) where gcd(B, C) = 1, find P, Q such that:
/// A/(B*C) = P/B + Q/C
///
/// This uses: 1 = s*B + t*C (Bézout identity)
/// Then: A = A*s*B + A*t*C
/// So: A/(B*C) = (A*t)/B + (A*s)/C
///
/// After reducing: P = (A*t) mod B, Q = (A*s) mod C
pub fn partial_fraction_coprime<K: Field>(
    a: &DensePoly<K>,
    b: &DensePoly<K>,
    c: &DensePoly<K>,
) -> Option<(DensePoly<K>, DensePoly<K>)> {
    let (gcd, s, t) = poly_extended_gcd(b, c);

    // Check that b and c are coprime
    if gcd.degree() != 0 {
        return None;
    }

    // Scale by gcd inverse
    let gcd_inv = gcd.coeff(0).inv()?;

    // P = (A*t) mod B
    let at = a.mul(&t).scale(&gcd_inv);
    let (_, p) = poly_div_rem(&at, b);

    // Q = (A*s) mod C
    let as_ = a.mul(&s).scale(&gcd_inv);
    let (_, q) = poly_div_rem(&as_, c);

    Some((p, q))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_poly::dense::DensePoly;
    use tertius_rings::rationals::Q;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_single_factor() {
        // 1/(x+1)^2 should decompose to A/(x+1) + B/(x+1)^2
        let num = poly(&[1]);
        let den = poly(&[1, 2, 1]); // (x+1)^2
        let rf = RationalFunction::new(num, den);

        let factor = poly(&[1, 1]); // x+1
        let decomp = partial_fraction_decomposition(&rf, &[(factor, 2)]);

        assert!(decomp.polynomial_part.is_zero());
        // Should have terms with powers 1 and 2
        assert!(!decomp.terms.is_empty());
    }

    #[test]
    fn test_polynomial_part_extraction() {
        // (x^2 + 1) / (x + 1) = (x - 1) + 2/(x+1)
        let num = poly(&[1, 0, 1]); // x^2 + 1
        let den = poly(&[1, 1]); // x + 1
        let rf = RationalFunction::new(num, den.clone());

        let decomp = partial_fraction_decomposition(&rf, &[(den, 1)]);

        // Polynomial part: x - 1
        assert_eq!(decomp.polynomial_part.degree(), 1);
    }

    #[test]
    fn test_coprime_simple() {
        // 1/(x(x+1)) = A/x + B/(x+1)
        // 1 = A(x+1) + Bx
        // x=0: 1 = A
        // x=-1: 1 = -B, so B = -1
        // Result: 1/x - 1/(x+1)

        let a = poly(&[1]); // 1
        let b = poly(&[0, 1]); // x
        let c = poly(&[1, 1]); // x + 1

        let result = partial_fraction_coprime(&a, &b, &c);
        assert!(result.is_some());

        let (p, r) = result.unwrap();

        // P should be 1 (coefficient for 1/x)
        assert_eq!(p.degree(), 0);
        assert_eq!(p.coeff(0), q(1));

        // Q should be -1 (coefficient for 1/(x+1))
        assert_eq!(r.degree(), 0);
        assert_eq!(r.coeff(0), q(-1));
    }

    #[test]
    fn test_coprime_factors_reconstruction() {
        // Test that decomposition can be reconstructed
        // 1 / ((x-1)(x+1)) = 1/(x^2-1)

        let num = poly(&[1]);
        let den = poly(&[-1, 0, 1]); // x^2 - 1
        let rf = RationalFunction::new(num, den);

        let f1 = poly(&[-1, 1]); // x - 1
        let f2 = poly(&[1, 1]); // x + 1

        let decomp = partial_fraction_decomposition(&rf, &[(f1, 1), (f2, 1)]);

        // Reconstruct and verify
        let reconstructed = decomp.to_rational_function();

        // They should be equal
        assert_eq!(reconstructed, rf);
    }
}
