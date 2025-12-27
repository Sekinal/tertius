//! The ring of integers Z.

use crate::traits::{CommutativeRing, EuclideanDomain, IntegralDomain, OrderedRing, Ring};
use tertius_integers::Integer;

/// The ring of integers.
///
/// This is a wrapper around `tertius_integers::Integer` that implements
/// the algebraic traits.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Z(pub Integer);

impl Z {
    /// Creates a new integer.
    #[must_use]
    pub fn new(value: i64) -> Self {
        Self(Integer::new(value))
    }

    /// Returns the inner Integer.
    #[must_use]
    pub fn into_inner(self) -> Integer {
        self.0
    }

    /// Returns a reference to the inner Integer.
    #[must_use]
    pub fn as_inner(&self) -> &Integer {
        &self.0
    }
}

impl Ring for Z {
    fn zero() -> Self {
        Self(Integer::new(0))
    }

    fn one() -> Self {
        Self(Integer::new(1))
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

impl CommutativeRing for Z {}
impl IntegralDomain for Z {}

impl EuclideanDomain for Z {
    fn div_rem(&self, other: &Self) -> (Self, Self) {
        let q = self.0.clone() / other.0.clone();
        let r = self.0.clone() % other.0.clone();
        (Self(q), Self(r))
    }

    fn gcd(&self, other: &Self) -> Self {
        Self(self.0.gcd(&other.0))
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self) {
        // Extended Euclidean algorithm
        let mut old_r = self.clone();
        let mut r = other.clone();
        let mut old_s = Self::one();
        let mut s = Self::zero();
        let mut old_t = Self::zero();
        let mut t = Self::one();

        while !r.is_zero() {
            let (q, rem) = old_r.div_rem(&r);
            old_r = r;
            r = rem;

            let new_s = old_s.clone() - q.clone() * s.clone();
            old_s = s;
            s = new_s;

            let new_t = old_t.clone() - q * t.clone();
            old_t = t;
            t = new_t;
        }

        (old_r, old_s, old_t)
    }
}

impl OrderedRing for Z {
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }

    fn signum(&self) -> i8 {
        self.0.signum()
    }
}

// Implement arithmetic operations
impl std::ops::Add for Z {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl std::ops::Sub for Z {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl std::ops::Mul for Z {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl std::ops::Neg for Z {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl From<i64> for Z {
    fn from(value: i64) -> Self {
        Self::new(value)
    }
}

impl From<Integer> for Z {
    fn from(value: Integer) -> Self {
        Self(value)
    }
}

impl std::fmt::Display for Z {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ring_laws() {
        let a = Z::new(10);
        let b = Z::new(3);

        // Test zero and one
        assert!(Z::zero().is_zero());
        assert!(Z::one().is_one());

        // Test arithmetic
        assert_eq!((a.clone() + b.clone()).0.to_i64(), Some(13));
        assert_eq!((a.clone() * b.clone()).0.to_i64(), Some(30));
    }

    #[test]
    fn test_euclidean_domain() {
        let a = Z::new(17);
        let b = Z::new(5);

        let (q, r) = a.div_rem(&b);
        assert_eq!(q.0.to_i64(), Some(3));
        assert_eq!(r.0.to_i64(), Some(2));
    }

    #[test]
    fn test_extended_gcd() {
        let a = Z::new(48);
        let b = Z::new(18);

        let (g, x, y) = a.extended_gcd(&b);
        assert_eq!(g.0.to_i64(), Some(6));

        // Verify: gcd = a*x + b*y
        let check = a * x + b * y;
        assert_eq!(check.0.to_i64(), Some(6));
    }
}
