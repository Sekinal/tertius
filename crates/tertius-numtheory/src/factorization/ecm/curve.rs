//! Montgomery curve representation for ECM.
//!
//! We use curves of the form: By² = x³ + Ax² + x
//! In projective coordinates: BY²Z = X³ + AX²Z + XZ²
//!
//! Montgomery curves allow efficient point doubling and differential addition,
//! which is all we need for ECM.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::point::ProjectivePoint;

/// A Montgomery curve modulo n.
///
/// Curve equation: By² = x³ + Ax² + x (mod n)
/// We only store A since B is not needed for x-coordinate operations.
#[derive(Debug, Clone)]
pub struct MontgomeryCurve {
    /// The A parameter (actually we store (A+2)/4 for efficiency).
    pub a24: Integer,
    /// The modulus.
    pub n: Integer,
}

impl MontgomeryCurve {
    /// Create a new Montgomery curve with parameter A.
    pub fn new(a: Integer, n: Integer) -> Self {
        // Store (A+2)/4 for efficient doubling
        let a_plus_2 = (a + Integer::from(2i64)) % n.clone();
        let four_inv = mod_inverse(&Integer::from(4i64), &n).unwrap_or(Integer::one());
        let a24 = (a_plus_2 * four_inv) % n.clone();
        Self { a24, n }
    }

    /// Point doubling: Q = 2P using Montgomery's formulas.
    ///
    /// Using projective coordinates (X : Z), the doubling formula is:
    /// X₂ = (X + Z)²(X - Z)²
    /// Z₂ = 4XZ((X - Z)² + ((A+2)/4)·4XZ)
    pub fn double(&self, p: &ProjectivePoint) -> Option<ProjectivePoint> {
        // u = (X + Z)²
        let sum = mod_add(&p.x, &p.z, &self.n);
        let u = mod_mul(&sum, &sum, &self.n);

        // v = (X - Z)²
        let diff = mod_sub(&p.x, &p.z, &self.n);
        let v = mod_mul(&diff, &diff, &self.n);

        // X₂ = u × v
        let x2 = mod_mul(&u, &v, &self.n);

        // w = u - v = 4XZ
        let w = mod_sub(&u, &v, &self.n);

        // Z₂ = w × (v + a24 × w)
        let a24_w = mod_mul(&self.a24, &w, &self.n);
        let v_plus_a24w = mod_add(&v, &a24_w, &self.n);
        let z2 = mod_mul(&w, &v_plus_a24w, &self.n);

        Some(ProjectivePoint::new(x2, z2))
    }

    /// Differential addition: Q = P + R given P - R.
    ///
    /// For Montgomery curves, we can compute P + R given P, R, and P - R.
    /// This is the Montgomery ladder primitive.
    pub fn add(&self, p: &ProjectivePoint, r: &ProjectivePoint, diff: &ProjectivePoint) -> Option<ProjectivePoint> {
        // u = (Xₚ - Zₚ)(Xᵣ + Zᵣ)
        let p_diff = mod_sub(&p.x, &p.z, &self.n);
        let r_sum = mod_add(&r.x, &r.z, &self.n);
        let u = mod_mul(&p_diff, &r_sum, &self.n);

        // v = (Xₚ + Zₚ)(Xᵣ - Zᵣ)
        let p_sum = mod_add(&p.x, &p.z, &self.n);
        let r_diff = mod_sub(&r.x, &r.z, &self.n);
        let v = mod_mul(&p_sum, &r_diff, &self.n);

        // X = Z_{P-R} × (u + v)²
        let sum = mod_add(&u, &v, &self.n);
        let sum_sq = mod_mul(&sum, &sum, &self.n);
        let x = mod_mul(&diff.z, &sum_sq, &self.n);

        // Z = X_{P-R} × (u - v)²
        let diff_uv = mod_sub(&u, &v, &self.n);
        let diff_sq = mod_mul(&diff_uv, &diff_uv, &self.n);
        let z = mod_mul(&diff.x, &diff_sq, &self.n);

        Some(ProjectivePoint::new(x, z))
    }

    /// Scalar multiplication using Montgomery ladder.
    ///
    /// Computes [k]P for scalar k.
    pub fn scalar_mul(&self, p: &ProjectivePoint, k: &Integer, n: &Integer) -> Option<ProjectivePoint> {
        if k.is_zero() {
            return Some(ProjectivePoint::identity());
        }
        if k.is_one() {
            return Some(p.clone());
        }

        // Montgomery ladder: constant-time and only needs x-coordinate
        let bits = integer_to_bits(k);

        // R0 = P, R1 = 2P
        let mut r0 = p.clone();
        let mut r1 = self.double(p)?;

        // Iterate from second-highest bit down to 0
        for i in 1..bits.len() {
            if bits[i] {
                // R0 = R0 + R1, R1 = 2R1
                r0 = self.add(&r0, &r1, p)?;
                r1 = self.double(&r1)?;
            } else {
                // R1 = R0 + R1, R0 = 2R0
                r1 = self.add(&r0, &r1, p)?;
                r0 = self.double(&r0)?;
            }
        }

        Some(r0)
    }
}

/// Convert integer to bits (MSB first).
fn integer_to_bits(n: &Integer) -> Vec<bool> {
    if n.is_zero() {
        return vec![false];
    }

    let mut bits = Vec::new();
    let mut m = n.clone();
    let two = Integer::from(2i64);

    while !m.is_zero() {
        bits.push(!(m.clone() % two.clone()).is_zero());
        m = m / two.clone();
    }

    bits.reverse();
    bits
}

/// Modular addition.
fn mod_add(a: &Integer, b: &Integer, n: &Integer) -> Integer {
    ((a.clone() % n.clone()) + (b.clone() % n.clone())) % n.clone()
}

/// Modular subtraction.
fn mod_sub(a: &Integer, b: &Integer, n: &Integer) -> Integer {
    let a = a.clone() % n.clone();
    let b = b.clone() % n.clone();
    if a >= b {
        a - b
    } else {
        a + n.clone() - b
    }
}

/// Modular multiplication.
fn mod_mul(a: &Integer, b: &Integer, n: &Integer) -> Integer {
    ((a.clone() % n.clone()) * (b.clone() % n.clone())) % n.clone()
}

/// Modular inverse.
fn mod_inverse(a: &Integer, n: &Integer) -> Option<Integer> {
    let (g, x, _) = extended_gcd(a, n);
    if g.is_one() {
        Some(((x % n.clone()) + n.clone()) % n.clone())
    } else {
        None
    }
}

/// Extended GCD.
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
    fn test_doubling() {
        // Test on a small modulus
        let n = Integer::from(143i64); // 11 × 13
        let curve = MontgomeryCurve::new(Integer::from(3i64), n.clone());
        let p = ProjectivePoint::new(Integer::from(5i64), Integer::one());

        let p2 = curve.double(&p).unwrap();
        // Just check that it doesn't crash and produces a valid point
        assert!(!p2.z.is_zero() || !p2.x.is_zero());
    }
}
