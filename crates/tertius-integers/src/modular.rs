//! Modular arithmetic.
//!
//! This module provides integers modulo a prime, essential for
//! FFT-based polynomial multiplication and Gröbner basis computation.

use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::Integer;

/// A modular integer with a compile-time prime modulus.
///
/// This is optimized for small primes that fit in a u64.
/// All operations are performed modulo P.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct ModInt<const P: u64>(u64);

impl<const P: u64> ModInt<P> {
    /// Creates a new modular integer.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value % P)
    }

    /// Creates a modular integer from a signed value.
    #[must_use]
    pub fn from_signed(value: i64) -> Self {
        if value >= 0 {
            Self::new(value as u64)
        } else {
            Self(P - ((-value) as u64 % P) % P)
        }
    }

    /// Returns the value as a u64.
    #[must_use]
    pub const fn value(self) -> u64 {
        self.0
    }

    /// Returns the modulus.
    #[must_use]
    pub const fn modulus() -> u64 {
        P
    }

    /// Computes the modular inverse using extended Euclidean algorithm.
    ///
    /// Returns `None` if the inverse doesn't exist (when gcd(self, P) != 1).
    #[must_use]
    pub fn inv(self) -> Option<Self> {
        if self.0 == 0 {
            return None;
        }

        // Extended Euclidean algorithm
        let mut t = 0i64;
        let mut new_t = 1i64;
        let mut r = P as i64;
        let mut new_r = self.0 as i64;

        while new_r != 0 {
            let quotient = r / new_r;
            (t, new_t) = (new_t, t - quotient * new_t);
            (r, new_r) = (new_r, r - quotient * new_r);
        }

        if r > 1 {
            return None; // Not coprime
        }

        Some(Self::from_signed(t))
    }

    /// Computes self^exp using binary exponentiation.
    #[must_use]
    pub fn pow(self, mut exp: u64) -> Self {
        let mut base = self;
        let mut result = Self::one();

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exp >>= 1;
        }

        result
    }
}

impl<const P: u64> Zero for ModInt<P> {
    fn zero() -> Self {
        Self(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl<const P: u64> One for ModInt<P> {
    fn one() -> Self {
        Self(1)
    }

    fn is_one(&self) -> bool {
        self.0 == 1
    }
}

impl<const P: u64> fmt::Debug for ModInt<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (mod {})", self.0, P)
    }
}

impl<const P: u64> fmt::Display for ModInt<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const P: u64> Add for ModInt<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self((self.0 + rhs.0) % P)
    }
}

impl<const P: u64> Sub for ModInt<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self((P + self.0 - rhs.0) % P)
    }
}

impl<const P: u64> Mul for ModInt<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // Use u128 to avoid overflow
        Self(((self.0 as u128 * rhs.0 as u128) % P as u128) as u64)
    }
}

impl<const P: u64> Div for ModInt<P> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv().expect("division by non-invertible element")
    }
}

impl<const P: u64> Neg for ModInt<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.0 == 0 {
            Self(0)
        } else {
            Self(P - self.0)
        }
    }
}

impl<const P: u64> From<u64> for ModInt<P> {
    fn from(value: u64) -> Self {
        Self::new(value)
    }
}

impl<const P: u64> From<i64> for ModInt<P> {
    fn from(value: i64) -> Self {
        Self::from_signed(value)
    }
}

/// Common NTT-friendly primes.
pub mod primes {
    /// 2^23 * 7 * 17 + 1 = 998244353 (common competitive programming prime)
    pub const P998244353: u64 = 998_244_353;

    /// 2^24 * 73 + 1 = 1224736769
    pub const P1224736769: u64 = 1_224_736_769;

    /// 2^25 * 59 + 1 = 1977716737
    pub const P1977716737: u64 = 1_977_716_737;

    /// Large prime for general use: 2^61 - 1 (Mersenne prime)
    pub const MERSENNE_61: u64 = (1 << 61) - 1;
}

/// Type alias for the common NTT prime.
pub type Mod998244353 = ModInt<{ primes::P998244353 }>;

/// A modular integer with a runtime-determined modulus.
///
/// Use this when the modulus is not known at compile time.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct ModIntDyn {
    value: Integer,
    modulus: Integer,
}

impl ModIntDyn {
    /// Creates a new modular integer.
    ///
    /// # Panics
    ///
    /// Panics if modulus is zero.
    #[must_use]
    pub fn new(value: Integer, modulus: Integer) -> Self {
        assert!(!modulus.is_zero(), "modulus cannot be zero");
        let value = if value.is_negative() {
            let r = value.clone() % modulus.clone();
            if r.is_zero() {
                r
            } else {
                r + modulus.clone()
            }
        } else {
            value % modulus.clone()
        };
        Self { value, modulus }
    }

    /// Returns the value.
    #[must_use]
    pub fn value(&self) -> &Integer {
        &self.value
    }

    /// Returns the modulus.
    #[must_use]
    pub fn modulus(&self) -> &Integer {
        &self.modulus
    }
}

impl fmt::Debug for ModIntDyn {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (mod {})", self.value, self.modulus)
    }
}

impl fmt::Display for ModIntDyn {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Mod7 = ModInt<7>;

    #[test]
    fn test_basic_ops() {
        let a = Mod7::new(5);
        let b = Mod7::new(4);

        assert_eq!((a + b).value(), 2); // 5 + 4 = 9 ≡ 2 (mod 7)
        assert_eq!((a - b).value(), 1); // 5 - 4 = 1
        assert_eq!((a * b).value(), 6); // 5 * 4 = 20 ≡ 6 (mod 7)
    }

    #[test]
    fn test_inverse() {
        // 3 * 5 = 15 ≡ 1 (mod 7), so inv(3) = 5
        assert_eq!(Mod7::new(3).inv(), Some(Mod7::new(5)));

        // 0 has no inverse
        assert_eq!(Mod7::new(0).inv(), None);
    }

    #[test]
    fn test_pow() {
        let a = Mod7::new(3);
        assert_eq!(a.pow(0).value(), 1);
        assert_eq!(a.pow(1).value(), 3);
        assert_eq!(a.pow(2).value(), 2); // 9 mod 7 = 2
        assert_eq!(a.pow(6).value(), 1); // Fermat's little theorem: a^(p-1) ≡ 1
    }

    #[test]
    fn test_negative() {
        let a = Mod7::from_signed(-3);
        assert_eq!(a.value(), 4); // -3 ≡ 4 (mod 7)
    }
}
