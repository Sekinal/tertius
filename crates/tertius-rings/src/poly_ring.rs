//! Polynomial rings R[x].
//!
//! This module provides the polynomial ring structure over any ring R.
//! For high-performance polynomial arithmetic, see `tertius-poly`.

use crate::traits::{CommutativeRing, IntegralDomain, Ring};

/// A polynomial over a ring R.
///
/// This is a generic implementation suitable for small polynomials.
/// For high-performance operations on large polynomials, use
/// `tertius_poly::DensePoly` instead.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Polynomial<R: Ring> {
    /// Coefficients in ascending degree order: [a_0, a_1, a_2, ...]
    /// Invariant: trailing zeros are removed (except for the zero polynomial).
    coeffs: Vec<R>,
}

impl<R: Ring> Polynomial<R> {
    /// Creates a new polynomial from coefficients.
    ///
    /// Coefficients are given in ascending degree order.
    #[must_use]
    pub fn new(mut coeffs: Vec<R>) -> Self {
        // Remove trailing zeros
        while coeffs.len() > 1 && coeffs.last().map_or(false, |c| c.is_zero()) {
            coeffs.pop();
        }

        if coeffs.is_empty() {
            coeffs.push(R::zero());
        }

        Self { coeffs }
    }

    /// Creates the zero polynomial.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            coeffs: vec![R::zero()],
        }
    }

    /// Creates the constant polynomial 1.
    #[must_use]
    pub fn one() -> Self {
        Self {
            coeffs: vec![R::one()],
        }
    }

    /// Creates a constant polynomial.
    #[must_use]
    pub fn constant(c: R) -> Self {
        Self::new(vec![c])
    }

    /// Creates the polynomial x.
    #[must_use]
    pub fn x() -> Self {
        Self::new(vec![R::zero(), R::one()])
    }

    /// Creates the monomial c * x^n.
    #[must_use]
    pub fn monomial(c: R, n: usize) -> Self {
        let mut coeffs = vec![R::zero(); n + 1];
        coeffs[n] = c;
        Self::new(coeffs)
    }

    /// Returns the degree of the polynomial.
    ///
    /// The zero polynomial has degree 0 by convention.
    #[must_use]
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            self.coeffs.len() - 1
        }
    }

    /// Returns true if this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0].is_zero()
    }

    /// Returns the leading coefficient.
    #[must_use]
    pub fn leading_coeff(&self) -> &R {
        self.coeffs.last().unwrap()
    }

    /// Returns the coefficient of x^i.
    #[must_use]
    pub fn coeff(&self, i: usize) -> R {
        self.coeffs.get(i).cloned().unwrap_or_else(R::zero)
    }

    /// Returns all coefficients.
    #[must_use]
    pub fn coeffs(&self) -> &[R] {
        &self.coeffs
    }

    /// Evaluates the polynomial at a point.
    #[must_use]
    pub fn eval(&self, x: &R) -> R {
        // Horner's method
        let mut result = R::zero();
        for c in self.coeffs.iter().rev() {
            result = result * x.clone() + c.clone();
        }
        result
    }

    /// Adds two polynomials.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut result = Vec::with_capacity(len);

        for i in 0..len {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            let b = other.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            result.push(a + b);
        }

        Self::new(result)
    }

    /// Negates a polynomial.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self::new(self.coeffs.iter().map(|c| -c.clone()).collect())
    }

    /// Subtracts two polynomials.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiplies two polynomials (schoolbook algorithm).
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }

        let n = self.coeffs.len();
        let m = other.coeffs.len();
        let mut result = vec![R::zero(); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                result[i + j] = result[i + j].clone() + self.coeffs[i].clone() * other.coeffs[j].clone();
            }
        }

        Self::new(result)
    }

    /// Multiplies by a scalar.
    #[must_use]
    pub fn scale(&self, c: &R) -> Self {
        Self::new(self.coeffs.iter().map(|x| x.clone() * c.clone()).collect())
    }

    /// Computes the formal derivative.
    #[must_use]
    pub fn derivative(&self) -> Self {
        if self.degree() == 0 {
            return Self::zero();
        }

        let mut result = Vec::with_capacity(self.coeffs.len() - 1);
        for (i, c) in self.coeffs.iter().skip(1).enumerate() {
            result.push(c.mul_by_scalar((i + 1) as i64));
        }

        Self::new(result)
    }
}

impl<R: Ring> std::fmt::Display for Polynomial<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut terms = Vec::new();
        for (i, c) in self.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }

            let term = match i {
                0 => format!("{c:?}"),
                1 => format!("{c:?}*x"),
                _ => format!("{c:?}*x^{i}"),
            };
            terms.push(term);
        }

        write!(f, "{}", terms.join(" + "))
    }
}

impl<R: CommutativeRing> Ring for Polynomial<R> {
    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::one()
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn is_one(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0].is_one()
    }
}

impl<R: CommutativeRing> std::ops::Add for Polynomial<R> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Polynomial::add(&self, &rhs)
    }
}

impl<R: CommutativeRing> std::ops::Sub for Polynomial<R> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Polynomial::sub(&self, &rhs)
    }
}

impl<R: CommutativeRing> std::ops::Mul for Polynomial<R> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Polynomial::mul(&self, &rhs)
    }
}

impl<R: CommutativeRing> std::ops::Neg for Polynomial<R> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Polynomial::neg(&self)
    }
}

impl<R: CommutativeRing> CommutativeRing for Polynomial<R> {}
impl<R: IntegralDomain> IntegralDomain for Polynomial<R> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rationals::Q;

    #[test]
    fn test_polynomial_basic() {
        // p(x) = 1 + 2x + 3x^2
        let p = Polynomial::new(vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)]);
        assert_eq!(p.degree(), 2);
        assert_eq!(p.coeff(0), Q::from_integer(1));
        assert_eq!(p.coeff(1), Q::from_integer(2));
        assert_eq!(p.coeff(2), Q::from_integer(3));
    }

    #[test]
    fn test_polynomial_arithmetic() {
        let p = Polynomial::new(vec![Q::from_integer(1), Q::from_integer(2)]); // 1 + 2x
        let q = Polynomial::new(vec![Q::from_integer(3), Q::from_integer(4)]); // 3 + 4x

        // (1 + 2x) + (3 + 4x) = 4 + 6x
        let sum = p.add(&q);
        assert_eq!(sum.coeff(0), Q::from_integer(4));
        assert_eq!(sum.coeff(1), Q::from_integer(6));

        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        let prod = p.mul(&q);
        assert_eq!(prod.coeff(0), Q::from_integer(3));
        assert_eq!(prod.coeff(1), Q::from_integer(10));
        assert_eq!(prod.coeff(2), Q::from_integer(8));
    }

    #[test]
    fn test_derivative() {
        // p(x) = 1 + 2x + 3x^2, p'(x) = 2 + 6x
        let p = Polynomial::new(vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)]);
        let dp = p.derivative();
        assert_eq!(dp.coeff(0), Q::from_integer(2));
        assert_eq!(dp.coeff(1), Q::from_integer(6));
    }
}
