//! Algebraic number fields.
//!
//! This module provides algebraic extensions of Q, represented
//! as Q(α) where α is a root of an irreducible polynomial.
//!
//! # Example
//!
//! ```
//! use tertius_rings::algebraic::{AlgebraicField, AlgebraicNumber};
//! use tertius_rings::rationals::Q;
//! use tertius_rings::traits::{Ring, Field};
//! use std::sync::Arc;
//!
//! // Create Q(√2)
//! let field = Arc::new(AlgebraicField::quadratic(2));
//! let sqrt2 = AlgebraicNumber::generator(Arc::clone(&field));
//!
//! // Compute (√2)^(-1) = √2/2
//! let inv = sqrt2.inv().unwrap();
//! let product = sqrt2.mul(&inv);
//! assert!(product.is_one());
//! ```

use std::ops::{Add, Mul, Neg, Sub};
use std::sync::Arc;

use crate::rationals::Q;
use crate::traits::{CommutativeRing, EuclideanDomain, Field, IntegralDomain, Ring};

/// An algebraic number field Q(α).
///
/// Elements are represented as polynomials of degree < n,
/// where n is the degree of the minimal polynomial.
#[derive(Clone, Debug)]
pub struct AlgebraicField {
    /// The minimal polynomial of the generator α.
    /// Stored as coefficients [a_0, a_1, ..., a_n] for a_0 + a_1*x + ... + a_n*x^n.
    min_poly: Vec<Q>,
    /// The degree of the extension.
    degree: usize,
    /// Pre-computed inverse of leading coefficient for fast reduction.
    lead_inv: Q,
}

impl AlgebraicField {
    /// Creates a new algebraic field from a minimal polynomial.
    ///
    /// The polynomial must be monic and irreducible.
    ///
    /// # Panics
    ///
    /// Panics if the polynomial has degree < 1.
    #[must_use]
    pub fn new(min_poly: Vec<Q>) -> Self {
        assert!(min_poly.len() >= 2, "minimal polynomial must have degree >= 1");

        let degree = min_poly.len() - 1;
        let lead_inv = min_poly.last().unwrap().inv().unwrap();

        Self {
            min_poly,
            degree,
            lead_inv,
        }
    }

    /// Creates the field Q(√d) for a square-free integer d.
    #[must_use]
    pub fn quadratic(d: i64) -> Self {
        // Minimal polynomial: x^2 - d
        Self::new(vec![Q::from_integer(-d), Q::zero(), Q::one()])
    }

    /// Returns the degree of the extension.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns the minimal polynomial.
    #[must_use]
    pub fn min_poly(&self) -> &[Q] {
        &self.min_poly
    }

    /// Reduces a polynomial modulo the minimal polynomial.
    ///
    /// Takes coefficients [a_0, a_1, ...] and reduces to degree < n.
    pub fn reduce(&self, coeffs: &mut Vec<Q>) {
        while coeffs.len() > self.degree {
            let high_coeff = coeffs.pop().unwrap();
            if high_coeff.is_zero() {
                continue;
            }

            // Subtract high_coeff * x^(len-1) ≡ -high_coeff * (min_poly - x^n) mod min_poly
            // This replaces x^n with -(a_0 + a_1*x + ... + a_{n-1}*x^{n-1})
            let scaled = high_coeff.clone() * self.lead_inv.clone();
            for (i, c) in self.min_poly.iter().take(self.degree).enumerate() {
                if coeffs.len() <= i {
                    coeffs.resize(i + 1, Q::zero());
                }
                coeffs[i] = coeffs[i].clone() - scaled.clone() * c.clone();
            }
        }

        // Remove trailing zeros
        while coeffs.len() > 1 && coeffs.last().map_or(false, |c| c.is_zero()) {
            coeffs.pop();
        }
    }
}

/// An element of an algebraic number field.
#[derive(Clone, Debug)]
pub struct AlgebraicNumber {
    /// Coefficients: a_0 + a_1*α + ... + a_{n-1}*α^{n-1}.
    coeffs: Vec<Q>,
    /// The field this element belongs to.
    field: Arc<AlgebraicField>,
}

impl AlgebraicNumber {
    /// Creates a new algebraic number.
    #[must_use]
    pub fn new(mut coeffs: Vec<Q>, field: Arc<AlgebraicField>) -> Self {
        field.reduce(&mut coeffs);
        Self { coeffs, field }
    }

    /// Creates the zero element.
    #[must_use]
    pub fn zero(field: Arc<AlgebraicField>) -> Self {
        Self {
            coeffs: vec![Q::zero()],
            field,
        }
    }

    /// Creates the one element.
    #[must_use]
    pub fn one(field: Arc<AlgebraicField>) -> Self {
        Self {
            coeffs: vec![Q::one()],
            field,
        }
    }

    /// Creates the generator α.
    #[must_use]
    pub fn generator(field: Arc<AlgebraicField>) -> Self {
        Self {
            coeffs: vec![Q::zero(), Q::one()],
            field,
        }
    }

    /// Creates an element from a rational.
    #[must_use]
    pub fn from_rational(r: Q, field: Arc<AlgebraicField>) -> Self {
        Self {
            coeffs: vec![r],
            field,
        }
    }

    /// Returns the coefficients.
    #[must_use]
    pub fn coeffs(&self) -> &[Q] {
        &self.coeffs
    }

    /// Returns the field.
    #[must_use]
    pub fn field(&self) -> &Arc<AlgebraicField> {
        &self.field
    }

    /// Returns true if this is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.is_zero())
    }

    /// Adds two algebraic numbers.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        assert!(Arc::ptr_eq(&self.field, &other.field), "fields must match");

        let len = self.coeffs.len().max(other.coeffs.len());
        let mut result = Vec::with_capacity(len);

        for i in 0..len {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(Q::zero);
            let b = other.coeffs.get(i).cloned().unwrap_or_else(Q::zero);
            result.push(a + b);
        }

        Self::new(result, Arc::clone(&self.field))
    }

    /// Negates an algebraic number.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(|c| -c.clone()).collect(),
            field: Arc::clone(&self.field),
        }
    }

    /// Subtracts two algebraic numbers.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiplies two algebraic numbers.
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        assert!(Arc::ptr_eq(&self.field, &other.field), "fields must match");

        let n = self.coeffs.len();
        let m = other.coeffs.len();
        let mut result = vec![Q::zero(); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                result[i + j] = result[i + j].clone() + self.coeffs[i].clone() * other.coeffs[j].clone();
            }
        }

        Self::new(result, Arc::clone(&self.field))
    }
}

impl PartialEq for AlgebraicNumber {
    fn eq(&self, other: &Self) -> bool {
        if !Arc::ptr_eq(&self.field, &other.field) {
            return false;
        }

        let max_len = self.coeffs.len().max(other.coeffs.len());
        for i in 0..max_len {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(Q::zero);
            let b = other.coeffs.get(i).cloned().unwrap_or_else(Q::zero);
            if a != b {
                return false;
            }
        }
        true
    }
}

impl Eq for AlgebraicNumber {}

impl std::fmt::Display for AlgebraicNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut terms = Vec::new();

        for (i, c) in self.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }

            let term = match i {
                0 => format!("{c}"),
                1 => {
                    if c.is_one() {
                        "α".to_string()
                    } else {
                        format!("{c}*α")
                    }
                }
                _ => {
                    if c.is_one() {
                        format!("α^{i}")
                    } else {
                        format!("{c}*α^{i}")
                    }
                }
            };
            terms.push(term);
        }

        if terms.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", terms.join(" + "))
        }
    }
}

// ============================================================================
// Polynomial helpers for computing inverses via extended GCD
// ============================================================================

/// Helper: Polynomial division with remainder over Q.
///
/// Given polynomials a(x) and b(x), computes (q, r) such that a = b*q + r
/// with deg(r) < deg(b).
fn poly_div_rem_q(a: &[Q], b: &[Q]) -> (Vec<Q>, Vec<Q>) {
    // Handle edge cases
    if b.is_empty() || (b.len() == 1 && b[0].is_zero()) {
        panic!("division by zero polynomial");
    }

    let b_deg = poly_degree(b);
    let a_deg = poly_degree(a);

    if a_deg < b_deg {
        return (vec![Q::zero()], a.to_vec());
    }

    let b_lead_inv = b[b_deg].inv().expect("leading coefficient must be invertible");

    let mut quotient = vec![Q::zero(); a_deg - b_deg + 1];
    let mut remainder = a.to_vec();

    while poly_degree(&remainder) >= b_deg && !poly_is_zero(&remainder) {
        let r_deg = poly_degree(&remainder);
        let deg_diff = r_deg - b_deg;
        let coeff = remainder[r_deg].clone() * b_lead_inv.clone();

        quotient[deg_diff] = coeff.clone();

        for (i, bc) in b.iter().enumerate().take(b_deg + 1) {
            if deg_diff + i < remainder.len() {
                remainder[deg_diff + i] = remainder[deg_diff + i].clone() - coeff.clone() * bc.clone();
            }
        }

        // Normalize
        poly_normalize(&mut remainder);
    }

    (quotient, remainder)
}

/// Helper: Extended GCD for polynomials over Q.
///
/// Returns (gcd, s, t) such that gcd = s*a + t*b.
/// The gcd is made monic.
fn poly_extended_gcd_q(a: &[Q], b: &[Q]) -> (Vec<Q>, Vec<Q>, Vec<Q>) {
    // Base cases
    if poly_is_zero(a) {
        if poly_is_zero(b) {
            return (vec![Q::zero()], vec![Q::one()], vec![Q::zero()]);
        }
        let g = poly_make_monic(b);
        let scale = b[poly_degree(b)].inv().expect("leading coeff should have inverse");
        return (g, vec![Q::zero()], vec![scale]);
    }
    if poly_is_zero(b) {
        let g = poly_make_monic(a);
        let scale = a[poly_degree(a)].inv().expect("leading coeff should have inverse");
        return (g, vec![scale], vec![Q::zero()]);
    }

    // Extended Euclidean algorithm
    let mut old_r = a.to_vec();
    let mut r = b.to_vec();
    let mut old_s = vec![Q::one()];
    let mut s = vec![Q::zero()];
    let mut old_t = vec![Q::zero()];
    let mut t = vec![Q::one()];

    while !poly_is_zero(&r) {
        let (q, rem) = poly_div_rem_q(&old_r, &r);

        let new_r = rem;
        let new_s = poly_sub_q(&old_s, &poly_mul_q(&q, &s));
        let new_t = poly_sub_q(&old_t, &poly_mul_q(&q, &t));

        old_r = r;
        r = new_r;
        old_s = s;
        s = new_s;
        old_t = t;
        t = new_t;
    }

    // Make the gcd monic and adjust s, t accordingly
    if poly_is_zero(&old_r) {
        return (vec![Q::zero()], vec![Q::one()], vec![Q::zero()]);
    }

    let lead = old_r[poly_degree(&old_r)].clone();
    let lead_inv = lead.inv().expect("leading coeff should have inverse");

    let gcd = poly_scale_q(&old_r, &lead_inv);
    let s_final = poly_scale_q(&old_s, &lead_inv);
    let t_final = poly_scale_q(&old_t, &lead_inv);

    (gcd, s_final, t_final)
}

/// Helper: Get polynomial degree.
fn poly_degree(p: &[Q]) -> usize {
    if p.is_empty() {
        return 0;
    }
    for i in (0..p.len()).rev() {
        if !p[i].is_zero() {
            return i;
        }
    }
    0
}

/// Helper: Check if polynomial is zero.
fn poly_is_zero(p: &[Q]) -> bool {
    p.is_empty() || p.iter().all(|c| c.is_zero())
}

/// Helper: Normalize polynomial (remove trailing zeros).
fn poly_normalize(p: &mut Vec<Q>) {
    while p.len() > 1 && p.last().map_or(false, |c| c.is_zero()) {
        p.pop();
    }
    if p.is_empty() {
        p.push(Q::zero());
    }
}

/// Helper: Make polynomial monic.
fn poly_make_monic(p: &[Q]) -> Vec<Q> {
    if poly_is_zero(p) {
        return vec![Q::zero()];
    }
    let deg = poly_degree(p);
    let lead_inv = p[deg].inv().expect("leading coeff should have inverse");
    poly_scale_q(p, &lead_inv)
}

/// Helper: Scale polynomial by scalar.
fn poly_scale_q(p: &[Q], c: &Q) -> Vec<Q> {
    p.iter().map(|x| x.clone() * c.clone()).collect()
}

/// Helper: Multiply two polynomials.
fn poly_mul_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    if poly_is_zero(a) || poly_is_zero(b) {
        return vec![Q::zero()];
    }

    let n = a.len();
    let m = b.len();
    let mut result = vec![Q::zero(); n + m - 1];

    for i in 0..n {
        for j in 0..m {
            result[i + j] = result[i + j].clone() + a[i].clone() * b[j].clone();
        }
    }

    poly_normalize(&mut result);
    result
}

/// Helper: Subtract two polynomials.
fn poly_sub_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let len = a.len().max(b.len());
    let mut result = Vec::with_capacity(len);

    for i in 0..len {
        let ai = a.get(i).cloned().unwrap_or_else(Q::zero);
        let bi = b.get(i).cloned().unwrap_or_else(Q::zero);
        result.push(ai - bi);
    }

    poly_normalize(&mut result);
    result
}

// ============================================================================
// Trait implementations for AlgebraicNumber
// ============================================================================

impl Add for AlgebraicNumber {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        AlgebraicNumber::add(&self, &rhs)
    }
}

impl Sub for AlgebraicNumber {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        AlgebraicNumber::sub(&self, &rhs)
    }
}

impl Mul for AlgebraicNumber {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        AlgebraicNumber::mul(&self, &rhs)
    }
}

impl Neg for AlgebraicNumber {
    type Output = Self;

    fn neg(self) -> Self::Output {
        AlgebraicNumber::neg(&self)
    }
}

impl Ring for AlgebraicNumber {
    fn zero() -> Self {
        panic!("AlgebraicNumber::zero() requires a field context; use AlgebraicNumber::zero(field) instead")
    }

    fn one() -> Self {
        panic!("AlgebraicNumber::one() requires a field context; use AlgebraicNumber::one(field) instead")
    }

    fn is_zero(&self) -> bool {
        AlgebraicNumber::is_zero(self)
    }

    fn is_one(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0].is_one()
    }
}

impl CommutativeRing for AlgebraicNumber {}
impl IntegralDomain for AlgebraicNumber {}

impl EuclideanDomain for AlgebraicNumber {
    fn div_rem(&self, other: &Self) -> (Self, Self) {
        // In a field, division is exact: a / b = a * b^(-1), remainder = 0
        let quotient = self.field_div(other);
        let zero = AlgebraicNumber::zero(Arc::clone(&self.field));
        (quotient, zero)
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self) {
        // In a field, gcd(a, b) = 1 for any non-zero elements
        // We can express 1 = a * (1/a) + b * 0 or similar
        if !self.is_zero() {
            let one = AlgebraicNumber::one(Arc::clone(&self.field));
            let zero = AlgebraicNumber::zero(Arc::clone(&self.field));
            let s = self.inv().unwrap();
            (one, s, zero)
        } else if !other.is_zero() {
            let one = AlgebraicNumber::one(Arc::clone(&self.field));
            let zero = AlgebraicNumber::zero(Arc::clone(&self.field));
            let t = other.inv().unwrap();
            (one, zero, t)
        } else {
            let zero = AlgebraicNumber::zero(Arc::clone(&self.field));
            let one = AlgebraicNumber::one(Arc::clone(&self.field));
            (zero.clone(), one, zero)
        }
    }
}

impl Field for AlgebraicNumber {
    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Use extended GCD to find inverse in Q(α) = Q[x]/(m(x))
        //
        // We have: a(α) is our element, m(α) = 0
        // Find s(x) such that s(x)·a(x) + t(x)·m(x) = 1 (mod m)
        // Then a(α)^(-1) = s(α)
        //
        // Since m(x) is irreducible and a(x) is not divisible by m(x),
        // gcd(a(x), m(x)) = 1, so such s exists.

        let (gcd, s, _t) = poly_extended_gcd_q(&self.coeffs, self.field.min_poly());

        // gcd should be 1 (constant) since m is irreducible
        if gcd.len() != 1 || !gcd[0].is_one() {
            // This shouldn't happen for irreducible minimal polynomials
            return None;
        }

        Some(AlgebraicNumber::new(s, Arc::clone(&self.field)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quadratic_field() {
        // Q(√2)
        let field = Arc::new(AlgebraicField::quadratic(2));
        assert_eq!(field.degree(), 2);

        let sqrt2 = AlgebraicNumber::generator(Arc::clone(&field));
        let two = AlgebraicNumber::from_rational(Q::from_integer(2), Arc::clone(&field));

        // √2 * √2 = 2
        let result = AlgebraicNumber::mul(&sqrt2, &sqrt2);
        assert_eq!(result, two);
    }

    #[test]
    fn test_arithmetic() {
        let field = Arc::new(AlgebraicField::quadratic(2));

        let a = AlgebraicNumber::new(
            vec![Q::from_integer(1), Q::from_integer(2)], // 1 + 2√2
            Arc::clone(&field),
        );
        let b = AlgebraicNumber::new(
            vec![Q::from_integer(3), Q::from_integer(-1)], // 3 - √2
            Arc::clone(&field),
        );

        // (1 + 2√2) + (3 - √2) = 4 + √2
        let sum = AlgebraicNumber::add(&a, &b);
        assert_eq!(sum.coeffs[0], Q::from_integer(4));
        assert_eq!(sum.coeffs[1], Q::from_integer(1));
    }

    #[test]
    fn test_inverse_sqrt2() {
        // Q(√2)
        let field = Arc::new(AlgebraicField::quadratic(2));
        let sqrt2 = AlgebraicNumber::generator(Arc::clone(&field));

        // (√2)^(-1) = √2/2
        let inv = sqrt2.inv().expect("√2 should have an inverse");

        // Verify: √2 * (√2)^(-1) = 1
        let product = AlgebraicNumber::mul(&sqrt2, &inv);
        assert!(product.is_one(), "√2 * (√2)^(-1) should be 1, got {:?}", product);

        // The inverse should be (1/2)√2
        assert_eq!(inv.coeffs.len(), 2);
        assert!(inv.coeffs[0].is_zero());
        assert_eq!(inv.coeffs[1], Q::new(1, 2));
    }

    #[test]
    fn test_inverse_rational() {
        // Q(√2)
        let field = Arc::new(AlgebraicField::quadratic(2));
        let three = AlgebraicNumber::from_rational(Q::from_integer(3), Arc::clone(&field));

        // 3^(-1) = 1/3
        let inv = three.inv().expect("3 should have an inverse");

        // Verify: 3 * (1/3) = 1
        let product = AlgebraicNumber::mul(&three, &inv);
        assert!(product.is_one());

        // The inverse should be 1/3
        assert_eq!(inv.coeffs.len(), 1);
        assert_eq!(inv.coeffs[0], Q::new(1, 3));
    }

    #[test]
    fn test_inverse_complex_element() {
        // Q(√2)
        let field = Arc::new(AlgebraicField::quadratic(2));

        // a = 1 + √2
        let a = AlgebraicNumber::new(
            vec![Q::from_integer(1), Q::from_integer(1)],
            Arc::clone(&field),
        );

        // a^(-1) should satisfy a * a^(-1) = 1
        let inv = a.inv().expect("1 + √2 should have an inverse");
        let product = AlgebraicNumber::mul(&a, &inv);
        assert!(product.is_one(), "a * a^(-1) should be 1, got {}", product);

        // Verify: (1 + √2)^(-1) = (1 - √2) / (1 - 2) = (1 - √2) / (-1) = -1 + √2
        assert_eq!(inv.coeffs.len(), 2);
        assert_eq!(inv.coeffs[0], Q::from_integer(-1));
        assert_eq!(inv.coeffs[1], Q::from_integer(1));
    }

    #[test]
    fn test_inverse_zero() {
        // Q(√2)
        let field = Arc::new(AlgebraicField::quadratic(2));
        let zero = AlgebraicNumber::zero(Arc::clone(&field));

        // 0 should not have an inverse
        assert!(zero.inv().is_none());
    }

    #[test]
    fn test_field_div() {
        // Q(√2)
        let field = Arc::new(AlgebraicField::quadratic(2));

        let a = AlgebraicNumber::new(
            vec![Q::from_integer(2), Q::from_integer(2)], // 2 + 2√2
            Arc::clone(&field),
        );
        let b = AlgebraicNumber::new(
            vec![Q::from_integer(1), Q::from_integer(1)], // 1 + √2
            Arc::clone(&field),
        );

        // (2 + 2√2) / (1 + √2) = 2(1 + √2) / (1 + √2) = 2
        let quotient = a.field_div(&b);
        let two = AlgebraicNumber::from_rational(Q::from_integer(2), Arc::clone(&field));
        assert_eq!(quotient, two);
    }

    #[test]
    fn test_cubic_field() {
        // Q(∛2): minimal polynomial x³ - 2
        let field = Arc::new(AlgebraicField::new(vec![
            Q::from_integer(-2),
            Q::zero(),
            Q::zero(),
            Q::one(),
        ]));
        assert_eq!(field.degree(), 3);

        let cbrt2 = AlgebraicNumber::generator(Arc::clone(&field));
        let two = AlgebraicNumber::from_rational(Q::from_integer(2), Arc::clone(&field));

        // (∛2)³ = 2
        let cubed = AlgebraicNumber::mul(&AlgebraicNumber::mul(&cbrt2, &cbrt2), &cbrt2);
        assert_eq!(cubed, two);

        // (∛2)^(-1) should satisfy ∛2 * (∛2)^(-1) = 1
        let inv = cbrt2.inv().expect("∛2 should have an inverse");
        let product = AlgebraicNumber::mul(&cbrt2, &inv);
        assert!(product.is_one(), "∛2 * (∛2)^(-1) should be 1, got {}", product);
    }

    #[test]
    fn test_cyclotomic_4() {
        // Q(i): minimal polynomial x² + 1
        let field = Arc::new(AlgebraicField::new(vec![
            Q::from_integer(1),
            Q::zero(),
            Q::one(),
        ]));
        assert_eq!(field.degree(), 2);

        let i = AlgebraicNumber::generator(Arc::clone(&field));
        let neg_one = AlgebraicNumber::from_rational(Q::from_integer(-1), Arc::clone(&field));

        // i² = -1
        let i_squared = AlgebraicNumber::mul(&i, &i);
        assert_eq!(i_squared, neg_one);

        // i^(-1) = -i
        let i_inv = i.inv().expect("i should have an inverse");
        let expected = AlgebraicNumber::neg(&i);
        assert_eq!(i_inv, expected);

        // Verify: i * (-i) = -i² = -(-1) = 1
        let product = AlgebraicNumber::mul(&i, &i_inv);
        assert!(product.is_one());
    }

    #[test]
    fn test_cyclotomic_8() {
        // Q(ζ₈) where ζ₈ = e^(2πi/8) = (1 + i)/√2
        // Minimal polynomial: x⁴ + 1
        let field = Arc::new(AlgebraicField::new(vec![
            Q::from_integer(1),
            Q::zero(),
            Q::zero(),
            Q::zero(),
            Q::one(),
        ]));
        assert_eq!(field.degree(), 4);

        let zeta = AlgebraicNumber::generator(Arc::clone(&field));
        let neg_one = AlgebraicNumber::from_rational(Q::from_integer(-1), Arc::clone(&field));

        // ζ⁴ = -1
        let zeta2 = AlgebraicNumber::mul(&zeta, &zeta);
        let zeta3 = AlgebraicNumber::mul(&zeta2, &zeta);
        let zeta4 = AlgebraicNumber::mul(&zeta3, &zeta);
        assert_eq!(zeta4, neg_one);

        // ζ^(-1) should satisfy ζ * ζ^(-1) = 1
        let inv = zeta.inv().expect("ζ₈ should have an inverse");
        let product = AlgebraicNumber::mul(&zeta, &inv);
        assert!(product.is_one(), "ζ₈ * ζ₈^(-1) should be 1, got {}", product);

        // ζ^(-1) = ζ⁷ = -ζ³ (since ζ⁸ = 1, ζ⁴ = -1)
        let expected = AlgebraicNumber::neg(&zeta3);
        assert_eq!(inv, expected);
    }

    #[test]
    fn test_poly_helpers() {
        // Test polynomial division
        let a = vec![Q::from_integer(-1), Q::zero(), Q::one()]; // x² - 1
        let b = vec![Q::from_integer(-1), Q::one()]; // x - 1

        let (q, r) = poly_div_rem_q(&a, &b);
        // (x² - 1) / (x - 1) = x + 1, remainder 0
        assert_eq!(q, vec![Q::from_integer(1), Q::from_integer(1)]);
        assert!(poly_is_zero(&r) || (r.len() == 1 && r[0].is_zero()));
    }

    #[test]
    fn test_extended_gcd() {
        // gcd(x² + 1, x - 1) over Q should be 1 (coprime)
        let a = vec![Q::from_integer(1), Q::zero(), Q::one()]; // x² + 1
        let b = vec![Q::from_integer(-1), Q::one()]; // x - 1

        let (gcd, s, t) = poly_extended_gcd_q(&a, &b);

        // gcd should be 1
        assert_eq!(gcd.len(), 1);
        assert!(gcd[0].is_one());

        // Verify Bezout identity: s*a + t*b = gcd
        let sa = poly_mul_q(&s, &a);
        let tb = poly_mul_q(&t, &b);
        let sum_len = sa.len().max(tb.len());
        let mut check = vec![Q::zero(); sum_len];
        for i in 0..sum_len {
            let sai = sa.get(i).cloned().unwrap_or_else(Q::zero);
            let tbi = tb.get(i).cloned().unwrap_or_else(Q::zero);
            check[i] = sai + tbi;
        }
        poly_normalize(&mut check);
        assert_eq!(check, gcd);
    }
}
