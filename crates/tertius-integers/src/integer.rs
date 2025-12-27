//! Arbitrary precision integers.
//!
//! This module provides a wrapper around `dashu::Integer` with
//! convenient operations for symbolic computation.

use dashu::base::{Abs, BitTest, Gcd, Signed as DashuSigned};
use dashu::integer::IBig;
use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

/// An arbitrary precision integer.
///
/// This type wraps `dashu::IBig` and provides the operations
/// needed for polynomial arithmetic and symbolic computation.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Integer(IBig);

impl Integer {
    /// Creates a new integer from an i64.
    #[must_use]
    pub fn new(value: i64) -> Self {
        Self(IBig::from(value))
    }

    /// Creates an integer from a string in the given base.
    ///
    /// # Errors
    ///
    /// Returns an error if the string is not a valid integer.
    pub fn from_str_radix(s: &str, radix: u32) -> Result<Self, dashu::base::error::ParseError> {
        IBig::from_str_radix(s, radix).map(Self)
    }

    /// Returns the absolute value.
    #[must_use]
    pub fn abs(&self) -> Self {
        Self(self.0.clone().abs())
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

    /// Returns true if this integer is negative.
    #[must_use]
    pub fn is_negative(&self) -> bool {
        DashuSigned::is_negative(&self.0)
    }

    /// Returns the number of bits needed to represent this integer.
    #[must_use]
    pub fn bit_len(&self) -> usize {
        self.0.bit_len()
    }

    /// Computes the greatest common divisor.
    #[must_use]
    pub fn gcd(&self, other: &Self) -> Self {
        Self(IBig::from(self.0.clone().gcd(other.0.clone())))
    }

    /// Computes the least common multiple.
    #[must_use]
    pub fn lcm(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let g = self.gcd(other);
        Self(&self.0 / &g.0 * &other.0).abs()
    }

    /// Returns the inner `dashu::IBig`.
    #[must_use]
    pub fn into_inner(self) -> IBig {
        self.0
    }

    /// Returns a reference to the inner `dashu::IBig`.
    #[must_use]
    pub fn as_inner(&self) -> &IBig {
        &self.0
    }

    /// Attempts to convert to an i64.
    ///
    /// Returns `None` if the value doesn't fit in an i64.
    #[must_use]
    pub fn to_i64(&self) -> Option<i64> {
        self.0.clone().try_into().ok()
    }

    /// Computes self^exp for non-negative exp.
    #[must_use]
    pub fn pow(&self, exp: u32) -> Self {
        Self(self.0.pow(exp as usize))
    }
}

impl Zero for Integer {
    fn zero() -> Self {
        Self(IBig::ZERO)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for Integer {
    fn one() -> Self {
        Self(IBig::ONE)
    }

    fn is_one(&self) -> bool {
        self.0 == IBig::ONE
    }
}

impl fmt::Debug for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Integer({})", self.0)
    }
}

impl fmt::Display for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// Arithmetic operations
impl Add for Integer {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Add<&Integer> for Integer {
    type Output = Self;

    fn add(self, rhs: &Integer) -> Self::Output {
        Self(self.0 + &rhs.0)
    }
}

impl Add for &Integer {
    type Output = Integer;

    fn add(self, rhs: Self) -> Self::Output {
        Integer(&self.0 + &rhs.0)
    }
}

impl Sub for Integer {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Sub<&Integer> for Integer {
    type Output = Self;

    fn sub(self, rhs: &Integer) -> Self::Output {
        Self(self.0 - &rhs.0)
    }
}

impl Sub for &Integer {
    type Output = Integer;

    fn sub(self, rhs: Self) -> Self::Output {
        Integer(&self.0 - &rhs.0)
    }
}

impl Mul for Integer {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul<&Integer> for Integer {
    type Output = Self;

    fn mul(self, rhs: &Integer) -> Self::Output {
        Self(self.0 * &rhs.0)
    }
}

impl Mul for &Integer {
    type Output = Integer;

    fn mul(self, rhs: Self) -> Self::Output {
        Integer(&self.0 * &rhs.0)
    }
}

impl Div for Integer {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0 / rhs.0)
    }
}

impl Div<&Integer> for Integer {
    type Output = Self;

    fn div(self, rhs: &Integer) -> Self::Output {
        Self(self.0 / &rhs.0)
    }
}

impl Rem for Integer {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self(self.0 % rhs.0)
    }
}

impl Neg for Integer {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl Neg for &Integer {
    type Output = Integer;

    fn neg(self) -> Self::Output {
        Integer(-&self.0)
    }
}

impl From<i64> for Integer {
    fn from(value: i64) -> Self {
        Self::new(value)
    }
}

impl From<i32> for Integer {
    fn from(value: i32) -> Self {
        Self::new(value as i64)
    }
}

impl From<u64> for Integer {
    fn from(value: u64) -> Self {
        Self(IBig::from(value))
    }
}

impl From<IBig> for Integer {
    fn from(value: IBig) -> Self {
        Self(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_ops() {
        let a = Integer::new(10);
        let b = Integer::new(3);

        assert_eq!((a.clone() + b.clone()).to_i64(), Some(13));
        assert_eq!((a.clone() - b.clone()).to_i64(), Some(7));
        assert_eq!((a.clone() * b.clone()).to_i64(), Some(30));
        assert_eq!((a.clone() / b.clone()).to_i64(), Some(3));
        assert_eq!((a % b).to_i64(), Some(1));
    }

    #[test]
    fn test_gcd() {
        let a = Integer::new(48);
        let b = Integer::new(18);
        assert_eq!(a.gcd(&b).to_i64(), Some(6));
    }

    #[test]
    fn test_large_numbers() {
        let a = Integer::from_str_radix("123456789012345678901234567890", 10).unwrap();
        let b = Integer::from_str_radix("987654321098765432109876543210", 10).unwrap();
        let sum = a + b;
        assert_eq!(sum.to_string(), "1111111110111111111011111111100");
    }
}
