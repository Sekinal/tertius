//! Sparse multivariate polynomials.
//!
//! This module provides sparse polynomial representation for
//! efficient handling of polynomials with few non-zero terms.

use tertius_rings::traits::Ring;

use crate::monomial::PackedMonomial;
use crate::ordering::MonomialOrder;

/// A sparse multivariate polynomial.
///
/// Terms are stored as (monomial, coefficient) pairs, sorted by
/// the monomial ordering.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SparsePoly<R: Ring> {
    /// Terms in sorted order (by monomial ordering).
    terms: Vec<(PackedMonomial, R)>,
    /// Number of variables.
    num_vars: usize,
    /// Monomial ordering used for sorting.
    order: MonomialOrder,
}

impl<R: Ring> SparsePoly<R> {
    /// Creates a new polynomial from terms.
    ///
    /// Terms are automatically sorted and combined.
    #[must_use]
    pub fn new(terms: Vec<(PackedMonomial, R)>, num_vars: usize, order: MonomialOrder) -> Self {
        let mut poly = Self {
            terms,
            num_vars,
            order,
        };
        poly.normalize();
        poly
    }

    /// Creates the zero polynomial.
    #[must_use]
    pub fn zero(num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: Vec::new(),
            num_vars,
            order,
        }
    }

    /// Creates the constant polynomial 1.
    #[must_use]
    pub fn one(num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: vec![(PackedMonomial::one(num_vars), R::one())],
            num_vars,
            order,
        }
    }

    /// Creates a constant polynomial.
    #[must_use]
    pub fn constant(c: R, num_vars: usize, order: MonomialOrder) -> Self {
        if c.is_zero() {
            Self::zero(num_vars, order)
        } else {
            Self {
                terms: vec![(PackedMonomial::one(num_vars), c)],
                num_vars,
                order,
            }
        }
    }

    /// Creates a single variable x_i.
    #[must_use]
    pub fn var(i: usize, num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: vec![(PackedMonomial::var(i, num_vars), R::one())],
            num_vars,
            order,
        }
    }

    /// Returns true if this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of terms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Returns true if there are no terms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of variables.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the monomial ordering.
    #[must_use]
    pub fn order(&self) -> MonomialOrder {
        self.order
    }

    /// Returns the terms.
    #[must_use]
    pub fn terms(&self) -> &[(PackedMonomial, R)] {
        &self.terms
    }

    /// Returns the leading monomial.
    #[must_use]
    pub fn leading_monomial(&self) -> Option<&PackedMonomial> {
        self.terms.first().map(|(m, _)| m)
    }

    /// Returns the leading coefficient.
    #[must_use]
    pub fn leading_coeff(&self) -> Option<&R> {
        self.terms.first().map(|(_, c)| c)
    }

    /// Returns the leading term (monomial, coefficient).
    #[must_use]
    pub fn leading_term(&self) -> Option<&(PackedMonomial, R)> {
        self.terms.first()
    }

    /// Sorts terms and combines like terms.
    fn normalize(&mut self) {
        // Sort by monomial order (descending for leading term first)
        self.terms.sort_by(|a, b| {
            self.order.compare(&b.0, &a.0, self.num_vars)
        });

        // Combine like terms
        let mut i = 0;
        while i < self.terms.len() {
            let mut j = i + 1;
            while j < self.terms.len() && self.terms[i].0 == self.terms[j].0 {
                let c = self.terms.remove(j).1;
                self.terms[i].1 = self.terms[i].1.clone() + c;
            }
            if self.terms[i].1.is_zero() {
                self.terms.remove(i);
            } else {
                i += 1;
            }
        }
    }

    /// Adds two polynomials.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.num_vars, other.num_vars);
        assert!(self.order == other.order);

        let mut terms = self.terms.clone();
        terms.extend(other.terms.clone());

        Self::new(terms, self.num_vars, self.order)
    }

    /// Negates a polynomial.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            terms: self.terms.iter().map(|(m, c)| (*m, -c.clone())).collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Subtracts two polynomials.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiplies two polynomials (schoolbook algorithm).
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.num_vars, other.num_vars);
        assert!(self.order == other.order);

        if self.is_zero() || other.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        let mut terms = Vec::with_capacity(self.len() * other.len());

        for (m1, c1) in &self.terms {
            for (m2, c2) in &other.terms {
                let m = m1.mul(m2);
                let c = c1.clone() * c2.clone();
                terms.push((m, c));
            }
        }

        Self::new(terms, self.num_vars, self.order)
    }

    /// Multiplies by a scalar.
    #[must_use]
    pub fn scale(&self, c: &R) -> Self {
        if c.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        Self {
            terms: self.terms.iter().map(|(m, x)| (*m, x.clone() * c.clone())).collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Multiplies by a monomial.
    #[must_use]
    pub fn mul_monomial(&self, m: &PackedMonomial, c: &R) -> Self {
        if c.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        Self {
            terms: self.terms.iter().map(|(m2, c2)| (m.mul(m2), c2.clone() * c.clone())).collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Computes the total degree.
    #[must_use]
    pub fn total_degree(&self) -> u32 {
        self.terms
            .iter()
            .map(|(m, _)| m.total_degree(self.num_vars))
            .max()
            .unwrap_or(0)
    }
}

impl<R: Ring> std::fmt::Display for SparsePoly<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let terms: Vec<_> = self
            .terms
            .iter()
            .map(|(m, c)| {
                let mon = m.to_string(self.num_vars);
                if mon == "1" {
                    format!("{c:?}")
                } else {
                    format!("{c:?}*{mon}")
                }
            })
            .collect();

        write!(f, "{}", terms.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_sparse_basic() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let y = SparsePoly::<Q>::var(1, 2, order);

        // x + y
        let sum = x.add(&y);
        assert_eq!(sum.len(), 2);
    }

    #[test]
    fn test_sparse_mul() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let one = SparsePoly::constant(Q::from_integer(1), 2, order);

        // (x + 1)^2 = x^2 + 2x + 1
        let xp1 = x.add(&one);
        let sq = xp1.mul(&xp1);
        assert_eq!(sq.len(), 3);
    }
}
