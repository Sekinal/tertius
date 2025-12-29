//! Projective point representation for ECM.
//!
//! For Montgomery curves, we only need the x-coordinate, represented
//! in projective form as (X : Z) where x = X/Z.

use num_traits::{One, Zero};
use tertius_integers::Integer;

/// A point on a Montgomery curve in projective coordinates.
///
/// Represents the point with x-coordinate X/Z.
/// The point at infinity is represented by Z = 0.
#[derive(Debug, Clone)]
pub struct ProjectivePoint {
    /// The X coordinate.
    pub x: Integer,
    /// The Z coordinate.
    pub z: Integer,
}

impl ProjectivePoint {
    /// Create a new projective point.
    pub fn new(x: Integer, z: Integer) -> Self {
        Self { x, z }
    }

    /// Create the identity point (point at infinity).
    pub fn identity() -> Self {
        Self {
            x: Integer::one(),
            z: Integer::zero(),
        }
    }

    /// Check if this is the identity point.
    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    /// Normalize the point to have Z = 1 (if possible).
    ///
    /// Returns None if Z is not invertible (which indicates a factor was found).
    pub fn normalize(&self, n: &Integer) -> Option<Self> {
        if self.z.is_zero() {
            return Some(Self::identity());
        }

        let z_inv = mod_inverse(&self.z, n)?;
        let x_norm = (self.x.clone() * z_inv) % n.clone();
        Some(Self::new(x_norm, Integer::one()))
    }
}

/// Modular inverse using extended Euclidean algorithm.
fn mod_inverse(a: &Integer, n: &Integer) -> Option<Integer> {
    if a.is_zero() {
        return None;
    }

    let (g, x, _) = extended_gcd(a, n);
    if g.is_one() {
        Some(((x % n.clone()) + n.clone()) % n.clone())
    } else {
        None
    }
}

/// Extended Euclidean algorithm.
fn extended_gcd(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    if b.is_zero() {
        return (a.clone(), Integer::one(), Integer::zero());
    }

    let (g, x1, y1) = extended_gcd(b, &(a.clone() % b.clone()));
    let x = y1.clone();
    let y = x1 - (a.clone() / b.clone()) * y1;
    (g, x, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity() {
        let id = ProjectivePoint::identity();
        assert!(id.is_identity());
    }

    #[test]
    fn test_normalize() {
        let p = ProjectivePoint::new(Integer::from(6i64), Integer::from(2i64));
        let n = Integer::from(7i64);
        let norm = p.normalize(&n).unwrap();
        assert_eq!(norm.x, Integer::from(3i64)); // 6/2 = 3
        assert_eq!(norm.z, Integer::from(1i64));
    }
}
