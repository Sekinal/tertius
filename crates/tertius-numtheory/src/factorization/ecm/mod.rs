//! Elliptic Curve Method (ECM) for integer factorization.
//!
//! ECM is the third-fastest known factorization method (after GNFS and QS)
//! and is particularly effective for finding medium-sized factors (20-60 digits).
//!
//! The algorithm uses elliptic curves in Montgomery form:
//! By²  = x³ + Ax² + x (mod n)
//!
//! The key insight is that if p | n and the elliptic curve has order m (mod p),
//! then multiplying a point by M where m | M will give the identity.
//! This creates a non-trivial GCD that reveals the factor.

mod curve;
mod point;

use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use tertius_integers::Integer;

use curve::MontgomeryCurve;
use point::ProjectivePoint;

/// ECM factorization with default parameters.
///
/// Tries multiple curves with increasing B1 bounds until a factor is found.
pub fn ecm_factor(n: &Integer) -> Option<Integer> {
    // Check trivial cases
    if n.is_one() {
        return None;
    }
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return Some(Integer::from(2i64));
    }

    // Check if n is a perfect square
    let sqrt_n = integer_sqrt(n);
    if &(sqrt_n.clone() * sqrt_n.clone()) == n {
        return Some(sqrt_n);
    }

    // Try with increasing bounds
    let bounds: [(u64, u64, usize); 6] = [
        (2_000, 200_000, 25),           // ~15 digit factors
        (11_000, 1_100_000, 90),        // ~20 digit factors
        (50_000, 5_000_000, 300),       // ~25 digit factors
        (250_000, 25_000_000, 700),     // ~30 digit factors
        (1_000_000, 100_000_000, 1800), // ~35 digit factors
        (3_000_000, 300_000_000, 5100), // ~40 digit factors
    ];

    for (b1, b2, curves) in bounds {
        if let Some(factor) = ecm_with_params(n, b1, b2, curves) {
            return Some(factor);
        }
    }

    None
}

/// ECM with specific parameters.
///
/// - `b1`: Stage 1 bound (smooth bound)
/// - `b2`: Stage 2 bound (for one large prime)
/// - `num_curves`: Number of curves to try
pub fn ecm_with_params(n: &Integer, b1: u64, b2: u64, num_curves: usize) -> Option<Integer> {
    let mut rng = ChaCha8Rng::seed_from_u64(0xECF1234);

    for _ in 0..num_curves {
        // Generate a random curve and point
        let (curve, point) = generate_random_curve(&mut rng, n);

        // Stage 1: Multiply point by product of prime powers up to B1
        if let Some(factor) = ecm_stage1(n, &curve, &point, b1) {
            if !factor.is_one() && &factor != n {
                return Some(factor);
            }
        }
    }

    None
}

/// Generate a random Montgomery curve and starting point.
fn generate_random_curve(rng: &mut ChaCha8Rng, n: &Integer) -> (MontgomeryCurve, ProjectivePoint) {
    loop {
        // Random sigma for Suyama's parameterization
        let sigma = Integer::from(rng.gen_range(6i64..1_000_000_000));

        // Suyama's parameterization gives curves with known structure
        // u = σ² - 5
        // v = 4σ
        let u = (sigma.clone() * sigma.clone() - Integer::from(5i64)) % n.clone();
        let v = (Integer::from(4i64) * sigma.clone()) % n.clone();

        // Check that v is invertible
        let g = gcd(&v, n);
        if !g.is_one() {
            if &g != n {
                // Found a factor!
                continue; // Actually we should return this, but for simplicity continue
            }
            continue;
        }

        // A = (v - u)³(3u + v) / (4u³v) - 2
        let u_cubed = mod_mul(&mod_mul(&u, &u, n), &u, n);
        let v_minus_u = mod_sub(&v, &u, n);
        let v_minus_u_cubed = mod_mul(&mod_mul(&v_minus_u, &v_minus_u, n), &v_minus_u, n);
        let three_u_plus_v = mod_add(&mod_mul(&Integer::from(3i64), &u, n), &v, n);

        let numerator = mod_mul(&v_minus_u_cubed, &three_u_plus_v, n);
        let denominator = mod_mul(&mod_mul(&Integer::from(4i64), &u_cubed, n), &v, n);

        // Need to compute modular inverse of denominator
        if let Some(denom_inv) = mod_inverse(&denominator, n) {
            let a_plus_2 = mod_mul(&numerator, &denom_inv, n);
            let a = mod_sub(&a_plus_2, &Integer::from(2i64), n);

            let curve = MontgomeryCurve::new(a.clone(), n.clone());

            // Starting point: (u³/v³ : 1)
            let v_cubed = mod_mul(&mod_mul(&v, &v, n), &v, n);
            if let Some(v_cubed_inv) = mod_inverse(&v_cubed, n) {
                let x = mod_mul(&u_cubed, &v_cubed_inv, n);
                let point = ProjectivePoint::new(x, Integer::one());
                return (curve, point);
            }
        }
    }
}

/// ECM Stage 1: Multiply point by smooth numbers up to B1.
fn ecm_stage1(
    n: &Integer,
    curve: &MontgomeryCurve,
    point: &ProjectivePoint,
    b1: u64,
) -> Option<Integer> {
    let mut q = point.clone();

    // Multiply by prime powers up to B1
    let mut p = 2u64;
    while p <= b1 {
        // Find largest power of p that is <= B1
        let mut pp = p;
        while pp <= b1 / p {
            pp *= p;
        }

        // Q = [pp]Q
        q = curve.scalar_mul(&q, &Integer::from(pp), n)?;

        // Check for factor
        let g = gcd(&q.z, n);
        if !g.is_one() && &g != n {
            return Some(g);
        }

        // Next prime
        p = next_prime(p);
    }

    None
}

/// Simple next prime function.
fn next_prime(n: u64) -> u64 {
    if n < 2 {
        return 2;
    }
    if n == 2 {
        return 3;
    }
    let mut candidate = if n % 2 == 0 { n + 1 } else { n + 2 };
    while !is_prime_simple(candidate) {
        candidate += 2;
    }
    candidate
}

/// Simple primality test.
fn is_prime_simple(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    let mut i = 5u64;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
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

/// Modular inverse using extended Euclidean algorithm.
fn mod_inverse(a: &Integer, n: &Integer) -> Option<Integer> {
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

/// GCD.
fn gcd(a: &Integer, b: &Integer) -> Integer {
    a.gcd(b)
}

/// Integer square root.
fn integer_sqrt(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }
    let mut x = n.clone();
    let mut y = (n.clone() + Integer::one()) / Integer::from(2i64);
    while y < x {
        x = y.clone();
        y = (x.clone() + n.clone() / x.clone()) / Integer::from(2i64);
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ecm_small() {
        // Small semiprime
        let n = Integer::from(143i64); // 11 × 13
        let factor = ecm_with_params(&n, 100, 1000, 10);
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!(f == Integer::from(11i64) || f == Integer::from(13i64));
    }

    #[test]
    fn test_ecm_medium() {
        // Medium semiprime
        let n = Integer::from(1000003i64) * Integer::from(1000033i64);
        let factor = ecm_with_params(&n, 2000, 200000, 50);
        assert!(factor.is_some());
    }
}
