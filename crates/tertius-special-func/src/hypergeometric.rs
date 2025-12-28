//! Hypergeometric functions.
//!
//! The generalized hypergeometric function is:
//!
//! ₚFᵧ(a₁,...,aₚ; b₁,...,bᵧ; z) = Σ_{n=0}^∞ [(a₁)ₙ⋯(aₚ)ₙ / ((b₁)ₙ⋯(bᵧ)ₙ)] zⁿ/n!
//!
//! where (a)ₙ = a(a+1)⋯(a+n-1) is the Pochhammer symbol.
//!
//! The most common is the Gauss hypergeometric function:
//!
//! ₂F₁(a, b; c; z) = Σ_{n=0}^∞ [(a)ₙ(b)ₙ / (c)ₙ] zⁿ/n!
//!
//! # Special Cases
//!
//! Many elementary and special functions are hypergeometric:
//! - (1-z)^(-a) = ₂F₁(a, b; b; z)
//! - ln(1+z)/z = ₂F₁(1, 1; 2; -z)
//! - arcsin(z)/z = ₂F₁(1/2, 1/2; 3/2; z²)
//! - Elliptic integrals can be expressed as ₂F₁

use tertius_rings::rationals::Q;
use tertius_rings::traits::Field;

/// The Gauss hypergeometric function ₂F₁(a, b; c; z).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hypergeometric2F1 {
    /// Parameter a.
    pub a: String,
    /// Parameter b.
    pub b: String,
    /// Parameter c (must not be zero or negative integer).
    pub c: String,
    /// Argument z.
    pub z: String,
}

impl Hypergeometric2F1 {
    /// Creates ₂F₁(a, b; c; z).
    pub fn new(a: &str, b: &str, c: &str, z: &str) -> Self {
        Self {
            a: a.to_string(),
            b: b.to_string(),
            c: c.to_string(),
            z: z.to_string(),
        }
    }

    /// Returns the derivative d/dz ₂F₁(a,b;c;z) = (ab/c) ₂F₁(a+1,b+1;c+1;z).
    pub fn derivative(&self) -> String {
        format!(
            "({0}*{1}/{2}) * ₂F₁({0}+1, {1}+1; {2}+1; {3})",
            self.a, self.b, self.c, self.z
        )
    }

    /// Pfaff transformation: ₂F₁(a,b;c;z) = (1-z)^(-a) ₂F₁(a, c-b; c; z/(z-1))
    pub fn pfaff_transformation(&self) -> String {
        format!(
            "(1-{3})^(-{0}) * ₂F₁({0}, {2}-{1}; {2}; {3}/({3}-1))",
            self.a, self.b, self.c, self.z
        )
    }

    /// Euler transformation: ₂F₁(a,b;c;z) = (1-z)^(c-a-b) ₂F₁(c-a, c-b; c; z)
    pub fn euler_transformation(&self) -> String {
        format!(
            "(1-{3})^({2}-{0}-{1}) * ₂F₁({2}-{0}, {2}-{1}; {2}; {3})",
            self.a, self.b, self.c, self.z
        )
    }
}

impl std::fmt::Display for Hypergeometric2F1 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "₂F₁({}, {}; {}; {})", self.a, self.b, self.c, self.z)
    }
}

/// Computes the Pochhammer symbol (a)_n = a(a+1)...(a+n-1).
pub fn pochhammer(a: Q, n: usize) -> Q {
    if n == 0 {
        return Q::from_integer(1);
    }

    let mut result = a.clone();
    let mut curr = a;

    for _ in 1..n {
        curr = curr + Q::from_integer(1);
        result = result * curr.clone();
    }

    result
}

/// Computes the hypergeometric series to n terms.
///
/// ₂F₁(a, b; c; z) ≈ Σ_{k=0}^{n-1} [(a)_k(b)_k / (c)_k] z^k / k!
pub fn hypergeometric_series(a: Q, b: Q, c: Q, z: Q, num_terms: usize) -> Option<Q> {
    // Check that c is not a non-positive integer
    // (simplified check for now)

    let mut result = Q::from_integer(0);
    let mut z_power = Q::from_integer(1);
    let mut factorial = Q::from_integer(1);

    for k in 0..num_terms {
        if k > 0 {
            z_power = z_power * z.clone();
            factorial = factorial * Q::from_integer(k as i64);
        }

        let a_k = pochhammer(a.clone(), k);
        let b_k = pochhammer(b.clone(), k);
        let c_k = pochhammer(c.clone(), k);

        if c_k == Q::from_integer(0) {
            return None; // c is a non-positive integer ≤ k
        }

        let term = (a_k * b_k * z_power.clone()) * (c_k * factorial.clone()).inv()?;
        result = result + term;
    }

    Some(result)
}

/// Special values and reductions of ₂F₁.
pub struct HypergeometricReductions;

impl HypergeometricReductions {
    /// (1-z)^(-a) = ₂F₁(a, b; b; z)
    pub fn power_function(a: &str) -> String {
        format!("₂F₁({0}, 1; 1; z) = (1-z)^(-{0})", a)
    }

    /// ln(1+z)/z = ₂F₁(1, 1; 2; -z)
    pub fn log_over_z() -> String {
        "₂F₁(1, 1; 2; -z) = ln(1+z)/z".to_string()
    }

    /// arcsin(z)/z = ₂F₁(1/2, 1/2; 3/2; z²)
    pub fn arcsin_over_z() -> String {
        "₂F₁(1/2, 1/2; 3/2; z²) = arcsin(z)/z".to_string()
    }

    /// arctan(z)/z = ₂F₁(1/2, 1; 3/2; -z²)
    pub fn arctan_over_z() -> String {
        "₂F₁(1/2, 1; 3/2; -z²) = arctan(z)/z".to_string()
    }

    /// Complete elliptic K(k) = (π/2) ₂F₁(1/2, 1/2; 1; k²)
    pub fn elliptic_k() -> String {
        "K(k) = (π/2) ₂F₁(1/2, 1/2; 1; k²)".to_string()
    }

    /// Complete elliptic E(k) = (π/2) ₂F₁(-1/2, 1/2; 1; k²)
    pub fn elliptic_e() -> String {
        "E(k) = (π/2) ₂F₁(-1/2, 1/2; 1; k²)".to_string()
    }
}

/// The confluent hypergeometric function ₁F₁(a; c; z) = M(a, c, z).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hypergeometric1F1 {
    /// Parameter a.
    pub a: String,
    /// Parameter c.
    pub c: String,
    /// Argument z.
    pub z: String,
}

impl Hypergeometric1F1 {
    /// Creates ₁F₁(a; c; z).
    pub fn new(a: &str, c: &str, z: &str) -> Self {
        Self {
            a: a.to_string(),
            c: c.to_string(),
            z: z.to_string(),
        }
    }
}

impl std::fmt::Display for Hypergeometric1F1 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "₁F₁({}; {}; {})", self.a, self.c, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hypergeometric_creation() {
        let h = Hypergeometric2F1::new("1/2", "1/2", "1", "k²");
        assert_eq!(format!("{}", h), "₂F₁(1/2, 1/2; 1; k²)");
    }

    #[test]
    fn test_pochhammer() {
        // (1)_3 = 1 * 2 * 3 = 6
        let result = pochhammer(Q::from_integer(1), 3);
        assert_eq!(result, Q::from_integer(6));

        // (a)_0 = 1
        let result = pochhammer(Q::from_integer(5), 0);
        assert_eq!(result, Q::from_integer(1));

        // (2)_2 = 2 * 3 = 6
        let result = pochhammer(Q::from_integer(2), 2);
        assert_eq!(result, Q::from_integer(6));
    }

    #[test]
    fn test_hypergeometric_series() {
        // ₂F₁(1, 1; 1; 0) = 1
        let result = hypergeometric_series(
            Q::from_integer(1),
            Q::from_integer(1),
            Q::from_integer(1),
            Q::from_integer(0),
            10,
        );
        assert_eq!(result, Some(Q::from_integer(1)));
    }

    #[test]
    fn test_derivative() {
        let h = Hypergeometric2F1::new("a", "b", "c", "z");
        let deriv = h.derivative();
        assert!(deriv.contains("₂F₁"));
    }

    #[test]
    fn test_reductions() {
        let log = HypergeometricReductions::log_over_z();
        assert!(log.contains("ln"));
    }
}
