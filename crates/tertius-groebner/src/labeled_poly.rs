//! Labeled polynomials with signatures for F5-style algorithms.

use crate::monomial::PackedMonomial;
use crate::signature::Signature;
use std::cmp::Ordering;

/// A polynomial with an associated F5 signature.
///
/// The signature tracks the polynomial's derivation from the original generators,
/// enabling detection of reductions to zero before they happen.
#[derive(Clone, Debug)]
pub struct LabeledPoly<R> {
    /// The polynomial as a list of (coefficient, monomial) pairs, sorted by monomial descending.
    terms: Vec<(R, PackedMonomial)>,
    /// The F5 signature.
    pub signature: Signature,
    /// Sugar degree (for selection heuristics).
    pub sugar: u32,
}

impl<R: Clone + PartialEq> LabeledPoly<R> {
    /// Creates a new labeled polynomial.
    pub fn new(terms: Vec<(R, PackedMonomial)>, signature: Signature, sugar: u32) -> Self {
        Self {
            terms,
            signature,
            sugar,
        }
    }

    /// Creates a labeled polynomial for a generator.
    pub fn from_generator(
        terms: Vec<(R, PackedMonomial)>,
        index: usize,
        num_vars: usize,
    ) -> Self {
        let sugar = terms
            .first()
            .map(|(_, m)| m.total_degree())
            .unwrap_or(0);

        Self {
            terms,
            signature: Signature::generator(index, num_vars),
            sugar,
        }
    }

    /// Returns the terms of the polynomial.
    pub fn terms(&self) -> &[(R, PackedMonomial)] {
        &self.terms
    }

    /// Returns the leading term (coefficient, monomial), if non-zero.
    pub fn leading_term(&self) -> Option<&(R, PackedMonomial)> {
        self.terms.first()
    }

    /// Returns the leading monomial, if non-zero.
    pub fn leading_monomial(&self) -> Option<&PackedMonomial> {
        self.terms.first().map(|(_, m)| m)
    }

    /// Returns the leading coefficient, if non-zero.
    pub fn leading_coeff(&self) -> Option<&R> {
        self.terms.first().map(|(c, _)| c)
    }

    /// Returns true if this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of terms.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// Returns the degree (total degree of leading monomial).
    pub fn degree(&self) -> u32 {
        self.terms.first().map(|(_, m)| m.total_degree()).unwrap_or(0)
    }

    /// Multiplies this polynomial by a monomial, updating signature accordingly.
    pub fn mul_monomial(&self, m: &PackedMonomial) -> Self
    where
        R: Clone,
    {
        let new_terms: Vec<_> = self
            .terms
            .iter()
            .map(|(c, mono)| (c.clone(), mono.mul(m)))
            .collect();

        Self {
            terms: new_terms,
            signature: self.signature.mul(m),
            sugar: self.sugar + m.total_degree(),
        }
    }

    /// Compares by signature using POT ordering.
    pub fn cmp_signature(&self, other: &Self) -> Ordering {
        self.signature.cmp_pot(&other.signature)
    }
}

/// A sparse polynomial without signature tracking (for intermediate computations).
#[derive(Clone, Debug)]
pub struct SparsePoly<R> {
    /// Terms sorted by monomial descending (grevlex).
    terms: Vec<(R, PackedMonomial)>,
}

impl<R: Clone + PartialEq + num_traits::Zero> SparsePoly<R> {
    /// Creates a new sparse polynomial from terms.
    ///
    /// Terms should be sorted by monomial descending.
    pub fn new(terms: Vec<(R, PackedMonomial)>) -> Self {
        // Filter out zero terms
        let terms: Vec<_> = terms.into_iter().filter(|(c, _)| !c.is_zero()).collect();
        Self { terms }
    }

    /// Creates the zero polynomial.
    pub fn zero() -> Self {
        Self { terms: vec![] }
    }

    /// Returns the terms.
    pub fn terms(&self) -> &[(R, PackedMonomial)] {
        &self.terms
    }

    /// Returns the leading monomial.
    pub fn leading_monomial(&self) -> Option<&PackedMonomial> {
        self.terms.first().map(|(_, m)| m)
    }

    /// Returns the leading coefficient.
    pub fn leading_coeff(&self) -> Option<&R> {
        self.terms.first().map(|(c, _)| c)
    }

    /// Returns true if zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of terms.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// Returns the total degree.
    pub fn degree(&self) -> u32 {
        self.terms.first().map(|(_, m)| m.total_degree()).unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_labeled_poly_leading() {
        let m1 = PackedMonomial::new(&[2, 1]); // x^2*y
        let m2 = PackedMonomial::new(&[1, 1]); // x*y
        let m3 = PackedMonomial::new(&[0, 1]); // y

        // Sort by grevlex: x^2*y > x*y > y
        let mut terms = vec![(1i64, m1.clone()), (2i64, m2.clone()), (-1i64, m3)];
        terms.sort_by(|a, b| b.1.cmp_grevlex(&a.1));

        let poly = LabeledPoly::new(terms, Signature::generator(0, 2), 3);

        assert_eq!(poly.leading_monomial(), Some(&m1));
        assert_eq!(poly.leading_coeff(), Some(&1i64));
        assert_eq!(poly.degree(), 3);
    }

    #[test]
    fn test_labeled_poly_mul_monomial() {
        let m = PackedMonomial::new(&[1, 0]); // x
        let terms = vec![(1i64, PackedMonomial::new(&[1, 1]))]; // x*y

        let poly = LabeledPoly::new(terms, Signature::generator(0, 2), 2);
        let product = poly.mul_monomial(&m);

        // x * (x*y) = x^2*y
        assert_eq!(product.leading_monomial().unwrap().exponent(0), 2);
        assert_eq!(product.leading_monomial().unwrap().exponent(1), 1);
        assert_eq!(product.sugar, 3);

        // Signature should also be multiplied by x
        assert_eq!(product.signature.monomial.exponent(0), 1);
    }
}
