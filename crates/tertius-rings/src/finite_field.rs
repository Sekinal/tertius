//! Finite fields Z_p and extensions.

use crate::traits::{CommutativeRing, EuclideanDomain, Field, IntegralDomain, Ring};
use tertius_integers::ModInt;

/// A finite field Z_p for prime p.
///
/// This wraps `ModInt<P>` and implements the algebraic traits.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct FiniteField<const P: u64>(pub ModInt<P>);

impl<const P: u64> FiniteField<P> {
    /// Creates a new field element.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(ModInt::new(value))
    }

    /// Creates a field element from a signed value.
    #[must_use]
    pub fn from_signed(value: i64) -> Self {
        Self(ModInt::from_signed(value))
    }

    /// Returns the value.
    #[must_use]
    pub const fn value(self) -> u64 {
        self.0.value()
    }

    /// Returns the characteristic (the prime p).
    #[must_use]
    pub const fn characteristic() -> u64 {
        P
    }
}

impl<const P: u64> Ring for FiniteField<P> {
    fn zero() -> Self {
        Self(ModInt::new(0))
    }

    fn one() -> Self {
        Self(ModInt::new(1))
    }

    fn is_zero(&self) -> bool {
        self.0.value() == 0
    }

    fn is_one(&self) -> bool {
        self.0.value() == 1
    }
}

impl<const P: u64> CommutativeRing for FiniteField<P> {}
impl<const P: u64> IntegralDomain for FiniteField<P> {}

impl<const P: u64> EuclideanDomain for FiniteField<P> {
    fn div_rem(&self, other: &Self) -> (Self, Self) {
        // In a field, division is exact
        (Self(self.0 / other.0), Self::zero())
    }

    fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() && other.is_zero() {
            Self::zero()
        } else {
            Self::one()
        }
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self) {
        if self.is_zero() && other.is_zero() {
            return (Self::zero(), Self::zero(), Self::zero());
        }

        if !self.is_zero() {
            (Self::one(), Self(self.0.inv().unwrap()), Self::zero())
        } else {
            (Self::one(), Self::zero(), Self(other.0.inv().unwrap()))
        }
    }
}

impl<const P: u64> Field for FiniteField<P> {
    fn inv(&self) -> Option<Self> {
        self.0.inv().map(Self)
    }
}

impl<const P: u64> std::ops::Add for FiniteField<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<const P: u64> std::ops::Sub for FiniteField<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<const P: u64> std::ops::Mul for FiniteField<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl<const P: u64> std::ops::Neg for FiniteField<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl<const P: u64> From<u64> for FiniteField<P> {
    fn from(value: u64) -> Self {
        Self::new(value)
    }
}

impl<const P: u64> From<i64> for FiniteField<P> {
    fn from(value: i64) -> Self {
        Self::from_signed(value)
    }
}

impl<const P: u64> std::fmt::Display for FiniteField<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Type alias for GF(2), the field with two elements.
pub type GF2 = FiniteField<2>;

/// Type alias for the common NTT prime field.
pub type GF998244353 = FiniteField<998_244_353>;

#[cfg(test)]
mod tests {
    use super::*;

    type F7 = FiniteField<7>;

    #[test]
    fn test_field_ops() {
        let a = F7::new(5);
        let b = F7::new(4);

        assert_eq!((a + b).value(), 2);
        assert_eq!((a - b).value(), 1);
        assert_eq!((a * b).value(), 6);
    }

    #[test]
    fn test_inverse() {
        let a = F7::new(3);
        let inv = a.inv().unwrap();
        assert_eq!((a * inv).value(), 1);
    }

    #[test]
    fn test_division() {
        let a = F7::new(5);
        let b = F7::new(3);
        let c = a.field_div(&b);

        // Verify: c * b = a
        assert_eq!((c * b).value(), a.value());
    }
}
