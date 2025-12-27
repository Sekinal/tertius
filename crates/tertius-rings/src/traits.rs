//! Algebraic structure traits.
//!
//! This module defines the core algebraic traits that form the foundation
//! of the type system for polynomials and other algebraic objects.

use std::fmt::Debug;
use std::ops::{Add, Mul, Neg, Sub};

/// A ring is a set with addition and multiplication operations.
///
/// # Laws
///
/// - Addition is associative and commutative with identity `zero()`
/// - Multiplication is associative with identity `one()`
/// - Multiplication distributes over addition
/// - Every element has an additive inverse (`neg`)
pub trait Ring:
    Clone + Eq + Debug + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
{
    /// The additive identity.
    fn zero() -> Self;

    /// The multiplicative identity.
    fn one() -> Self;

    /// Returns true if this is the additive identity.
    fn is_zero(&self) -> bool;

    /// Returns true if this is the multiplicative identity.
    fn is_one(&self) -> bool;

    /// Computes self + self + ... (n times).
    fn mul_by_scalar(&self, n: i64) -> Self {
        if n == 0 {
            return Self::zero();
        }

        let mut result = self.clone();
        let abs_n = n.unsigned_abs();

        for _ in 1..abs_n {
            result = result + self.clone();
        }

        if n < 0 {
            -result
        } else {
            result
        }
    }

    /// Computes self^n for non-negative n.
    fn pow(&self, n: u32) -> Self {
        if n == 0 {
            return Self::one();
        }

        let mut result = Self::one();
        let mut base = self.clone();
        let mut exp = n;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base.clone();
            }
            base = base.clone() * base;
            exp >>= 1;
        }

        result
    }
}

/// A commutative ring where multiplication is commutative.
///
/// Most rings in symbolic computation are commutative.
pub trait CommutativeRing: Ring {}

/// An integral domain is a commutative ring with no zero divisors.
///
/// If a * b = 0, then a = 0 or b = 0.
pub trait IntegralDomain: CommutativeRing {}

/// A Euclidean domain supports division with remainder.
///
/// For any a, b with b ≠ 0, there exist q, r such that:
/// - a = b*q + r
/// - Either r = 0 or φ(r) < φ(b) for some Euclidean function φ
pub trait EuclideanDomain: IntegralDomain {
    /// Computes the quotient and remainder of division.
    ///
    /// # Panics
    ///
    /// May panic if `other` is zero.
    fn div_rem(&self, other: &Self) -> (Self, Self);

    /// Computes the quotient of division.
    fn div(&self, other: &Self) -> Self {
        self.div_rem(other).0
    }

    /// Computes the remainder of division.
    fn rem(&self, other: &Self) -> Self {
        self.div_rem(other).1
    }

    /// Computes the greatest common divisor.
    fn gcd(&self, other: &Self) -> Self {
        let mut a = self.clone();
        let mut b = other.clone();

        while !b.is_zero() {
            let r = a.rem(&b);
            a = b;
            b = r;
        }

        a
    }

    /// Computes the least common multiple.
    fn lcm(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let g = self.gcd(other);
        self.div(&g) * other.clone()
    }

    /// Extended Euclidean algorithm.
    ///
    /// Returns (gcd, x, y) such that gcd = self*x + other*y.
    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self);
}

/// A field is a ring where every non-zero element has a multiplicative inverse.
pub trait Field: EuclideanDomain {
    /// Computes the multiplicative inverse.
    ///
    /// Returns `None` if the element is zero.
    fn inv(&self) -> Option<Self>;

    /// Divides by another element.
    ///
    /// # Panics
    ///
    /// Panics if `other` is zero.
    fn field_div(&self, other: &Self) -> Self {
        self.clone() * other.inv().expect("division by zero")
    }
}

/// Marker trait for ordered rings.
pub trait OrderedRing: Ring + Ord {
    /// Returns the absolute value.
    fn abs(&self) -> Self;

    /// Returns the sign: -1, 0, or 1.
    fn signum(&self) -> i8;
}

/// A unique factorization domain.
///
/// Every non-zero, non-unit element can be written as a product of
/// irreducible elements, unique up to order and associates.
pub trait UniqueFactorizationDomain: IntegralDomain {
    /// The type of irreducible factors.
    type Factor: Clone + Eq;

    /// Factorizes into irreducible factors.
    fn factor(&self) -> Vec<(Self::Factor, u32)>;

    /// Returns true if this element is irreducible.
    fn is_irreducible(&self) -> bool;
}

#[cfg(test)]
mod tests {
    // Tests will be added when concrete implementations are available
}
