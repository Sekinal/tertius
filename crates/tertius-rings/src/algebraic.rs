//! Algebraic number fields.
//!
//! This module provides algebraic extensions of Q, represented
//! as Q(α) where α is a root of an irreducible polynomial.

use std::sync::Arc;

use crate::rationals::Q;
use crate::traits::{Field, Ring};

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
        let result = sqrt2.mul(&sqrt2);
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
        let sum = a.add(&b);
        assert_eq!(sum.coeffs[0], Q::from_integer(4));
        assert_eq!(sum.coeffs[1], Q::from_integer(1));
    }
}
