//! Arbitrary precision rational numbers.
//!
//! This module provides exact rational arithmetic for symbolic computation.

use dashu::base::{Abs, Inverse, Signed as DashuSigned, UnsignedAbs};
use dashu::rational::RBig;
use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::Integer;

/// An arbitrary precision rational number.
///
/// Rationals are always stored in lowest terms with a positive denominator.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Rational(RBig);

impl Rational {
    /// Creates a new rational from numerator and denominator.
    ///
    /// # Panics
    ///
    /// Panics if the denominator is zero.
    #[must_use]
    pub fn new(numerator: Integer, denominator: Integer) -> Self {
        assert!(!denominator.is_zero(), "denominator cannot be zero");
        Self(RBig::from_parts(
            numerator.into_inner(),
            denominator.into_inner().unsigned_abs(),
        ))
    }

    /// Creates a rational from an integer (denominator = 1).
    #[must_use]
    pub fn from_integer(n: Integer) -> Self {
        Self(RBig::from(n.into_inner()))
    }

    /// Creates a rational from i64 numerator and denominator.
    ///
    /// # Panics
    ///
    /// Panics if the denominator is zero.
    #[must_use]
    pub fn from_i64(numerator: i64, denominator: i64) -> Self {
        Self::new(Integer::new(numerator), Integer::new(denominator))
    }

    /// Returns the numerator.
    #[must_use]
    pub fn numerator(&self) -> Integer {
        Integer::from(self.0.numerator().clone())
    }

    /// Returns the denominator.
    #[must_use]
    pub fn denominator(&self) -> Integer {
        Integer::from(dashu::integer::IBig::from(self.0.denominator().clone()))
    }

    /// Returns true if this rational is an integer.
    #[must_use]
    pub fn is_integer(&self) -> bool {
        self.0.denominator().is_one()
    }

    /// Converts to an integer if the denominator is 1.
    #[must_use]
    pub fn to_integer(&self) -> Option<Integer> {
        if self.is_integer() {
            Some(self.numerator())
        } else {
            None
        }
    }

    /// Returns the absolute value.
    #[must_use]
    pub fn abs(&self) -> Self {
        Self(self.0.clone().abs())
    }

    /// Returns the reciprocal (1/x).
    ///
    /// # Panics
    ///
    /// Panics if the rational is zero.
    #[must_use]
    pub fn recip(&self) -> Self {
        assert!(!self.is_zero(), "cannot take reciprocal of zero");
        Self(self.0.clone().inv())
    }

    /// Returns the sign: -1, 0, or 1.
    #[must_use]
    pub fn signum(&self) -> i8 {
        if self.0.is_zero() {
            0
        } else if DashuSigned::is_positive(&self.0) {
            1
        } else {
            -1
        }
    }

    /// Returns true if negative.
    #[must_use]
    pub fn is_negative(&self) -> bool {
        DashuSigned::is_negative(&self.0)
    }

    /// Returns the inner `dashu::RBig`.
    #[must_use]
    pub fn into_inner(self) -> RBig {
        self.0
    }

    /// Returns a reference to the inner `dashu::RBig`.
    #[must_use]
    pub fn as_inner(&self) -> &RBig {
        &self.0
    }

    /// Computes self^exp for non-negative exp.
    #[must_use]
    pub fn pow(&self, exp: u32) -> Self {
        Self(self.0.pow(exp as usize))
    }
}

impl Zero for Rational {
    fn zero() -> Self {
        Self(RBig::ZERO)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for Rational {
    fn one() -> Self {
        Self(RBig::ONE)
    }

    fn is_one(&self) -> bool {
        self.0 == RBig::ONE
    }
}

impl fmt::Debug for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Rational({})", self.0)
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_integer() {
            write!(f, "{}", self.numerator())
        } else {
            write!(f, "{}/{}", self.numerator(), self.denominator())
        }
    }
}

// Arithmetic operations
impl Add for Rational {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Add<&Rational> for Rational {
    type Output = Self;

    fn add(self, rhs: &Rational) -> Self::Output {
        Self(self.0 + &rhs.0)
    }
}

impl Add for &Rational {
    type Output = Rational;

    fn add(self, rhs: Self) -> Self::Output {
        Rational(&self.0 + &rhs.0)
    }
}

impl Sub for Rational {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Sub<&Rational> for Rational {
    type Output = Self;

    fn sub(self, rhs: &Rational) -> Self::Output {
        Self(self.0 - &rhs.0)
    }
}

impl Sub for &Rational {
    type Output = Rational;

    fn sub(self, rhs: Self) -> Self::Output {
        Rational(&self.0 - &rhs.0)
    }
}

impl Mul for Rational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul<&Rational> for Rational {
    type Output = Self;

    fn mul(self, rhs: &Rational) -> Self::Output {
        Self(self.0 * &rhs.0)
    }
}

impl Mul for &Rational {
    type Output = Rational;

    fn mul(self, rhs: Self) -> Self::Output {
        Rational(&self.0 * &rhs.0)
    }
}

impl Div for Rational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0 / rhs.0)
    }
}

impl Div<&Rational> for Rational {
    type Output = Self;

    fn div(self, rhs: &Rational) -> Self::Output {
        Self(self.0 / &rhs.0)
    }
}

impl Neg for Rational {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl Neg for &Rational {
    type Output = Rational;

    fn neg(self) -> Self::Output {
        Rational(-&self.0)
    }
}

impl From<Integer> for Rational {
    fn from(n: Integer) -> Self {
        Self::from_integer(n)
    }
}

impl From<i64> for Rational {
    fn from(n: i64) -> Self {
        Self::from_integer(Integer::new(n))
    }
}

impl From<i32> for Rational {
    fn from(n: i32) -> Self {
        Self::from_integer(Integer::new(n as i64))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_ops() {
        let a = Rational::from_i64(1, 2);
        let b = Rational::from_i64(1, 3);

        // 1/2 + 1/3 = 5/6
        let sum = a.clone() + b.clone();
        assert_eq!(sum.numerator().to_i64(), Some(5));
        assert_eq!(sum.denominator().to_i64(), Some(6));

        // 1/2 * 1/3 = 1/6
        let prod = a.clone() * b.clone();
        assert_eq!(prod.numerator().to_i64(), Some(1));
        assert_eq!(prod.denominator().to_i64(), Some(6));
    }

    #[test]
    fn test_reduction() {
        // 4/6 should reduce to 2/3
        let r = Rational::from_i64(4, 6);
        assert_eq!(r.numerator().to_i64(), Some(2));
        assert_eq!(r.denominator().to_i64(), Some(3));
    }

    #[test]
    fn test_display() {
        assert_eq!(Rational::from_i64(3, 1).to_string(), "3");
        assert_eq!(Rational::from_i64(2, 3).to_string(), "2/3");
    }
}
