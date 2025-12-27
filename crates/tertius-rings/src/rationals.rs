//! The field of rational numbers Q.

use crate::traits::{CommutativeRing, EuclideanDomain, Field, IntegralDomain, OrderedRing, Ring};
use tertius_integers::Rational;

/// The field of rational numbers.
///
/// This is a wrapper around `tertius_integers::Rational` that implements
/// the algebraic traits.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Q(pub Rational);

impl Q {
    /// Creates a new rational from numerator and denominator.
    #[must_use]
    pub fn new(num: i64, den: i64) -> Self {
        Self(Rational::from_i64(num, den))
    }

    /// Creates a rational from an integer.
    #[must_use]
    pub fn from_integer(n: i64) -> Self {
        Self(Rational::from(n))
    }

    /// Returns the inner Rational.
    #[must_use]
    pub fn into_inner(self) -> Rational {
        self.0
    }

    /// Returns a reference to the inner Rational.
    #[must_use]
    pub fn as_inner(&self) -> &Rational {
        &self.0
    }
}

impl Ring for Q {
    fn zero() -> Self {
        Self(Rational::from(0))
    }

    fn one() -> Self {
        Self(Rational::from(1))
    }

    fn is_zero(&self) -> bool {
        use num_traits::Zero;
        self.0.is_zero()
    }

    fn is_one(&self) -> bool {
        use num_traits::One;
        self.0.is_one()
    }
}

impl CommutativeRing for Q {}
impl IntegralDomain for Q {}

impl EuclideanDomain for Q {
    fn div_rem(&self, other: &Self) -> (Self, Self) {
        // In a field, division is exact, so remainder is always zero
        (Self(self.0.clone() / other.0.clone()), Self::zero())
    }

    fn gcd(&self, other: &Self) -> Self {
        // In a field, gcd of any two non-zero elements is 1
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
            // gcd = 1 = self * (1/self) + other * 0
            (
                Self::one(),
                Self(self.0.recip()),
                Self::zero(),
            )
        } else {
            // gcd = 1 = self * 0 + other * (1/other)
            (
                Self::one(),
                Self::zero(),
                Self(other.0.recip()),
            )
        }
    }
}

impl Field for Q {
    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Self(self.0.recip()))
        }
    }
}

impl OrderedRing for Q {
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }

    fn signum(&self) -> i8 {
        self.0.signum()
    }
}

// Implement arithmetic operations
impl std::ops::Add for Q {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl std::ops::Sub for Q {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl std::ops::Mul for Q {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl std::ops::Neg for Q {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl From<i64> for Q {
    fn from(value: i64) -> Self {
        Self::from_integer(value)
    }
}

impl From<Rational> for Q {
    fn from(value: Rational) -> Self {
        Self(value)
    }
}

impl std::fmt::Display for Q {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_laws() {
        let a = Q::new(2, 3);
        let b = Q::new(3, 4);

        // 2/3 + 3/4 = 8/12 + 9/12 = 17/12
        let sum = a.clone() + b.clone();
        assert_eq!(sum.0.numerator().to_i64(), Some(17));
        assert_eq!(sum.0.denominator().to_i64(), Some(12));

        // 2/3 * 3/4 = 6/12 = 1/2
        let prod = a.clone() * b.clone();
        assert_eq!(prod.0.numerator().to_i64(), Some(1));
        assert_eq!(prod.0.denominator().to_i64(), Some(2));
    }

    #[test]
    fn test_inverse() {
        let a = Q::new(3, 5);
        let inv = a.clone().inv().unwrap();

        // 3/5 * 5/3 = 1
        let prod = a * inv;
        assert!(prod.is_one());
    }

    #[test]
    fn test_division() {
        let a = Q::new(1, 2);
        let b = Q::new(1, 3);

        // (1/2) / (1/3) = (1/2) * 3 = 3/2
        let quot = a.field_div(&b);
        assert_eq!(quot.0.numerator().to_i64(), Some(3));
        assert_eq!(quot.0.denominator().to_i64(), Some(2));
    }
}
