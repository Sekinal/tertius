//! Elliptic integrals.
//!
//! The three standard elliptic integrals are:
//!
//! **Incomplete elliptic integral of the first kind:**
//! F(φ, k) = ∫₀^φ dθ / √(1 - k²sin²θ)
//!
//! **Incomplete elliptic integral of the second kind:**
//! E(φ, k) = ∫₀^φ √(1 - k²sin²θ) dθ
//!
//! **Incomplete elliptic integral of the third kind:**
//! Π(n; φ, k) = ∫₀^φ dθ / ((1 - n·sin²θ)√(1 - k²sin²θ))
//!
//! Complete elliptic integrals are obtained when φ = π/2:
//! K(k) = F(π/2, k), E(k) = E(π/2, k)
//!
//! # Integration
//!
//! Elliptic integrals arise from:
//! - ∫ 1/√(1-x⁴) dx → elliptic integral
//! - ∫ √(1-x²·k²) dx → E-type
//! - Arc length of ellipse

/// Incomplete elliptic integral of the first kind F(φ, k).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EllipticF {
    /// Amplitude φ (as symbolic expression).
    pub amplitude: String,
    /// Modulus k (as symbolic expression).
    pub modulus: String,
    /// Whether this is complete (φ = π/2).
    pub is_complete: bool,
}

impl EllipticF {
    /// Creates F(φ, k).
    pub fn incomplete(amplitude: &str, modulus: &str) -> Self {
        Self {
            amplitude: amplitude.to_string(),
            modulus: modulus.to_string(),
            is_complete: false,
        }
    }

    /// Creates the complete elliptic integral K(k) = F(π/2, k).
    pub fn complete(modulus: &str) -> Self {
        Self {
            amplitude: "π/2".to_string(),
            modulus: modulus.to_string(),
            is_complete: true,
        }
    }

    /// Returns the derivative with respect to amplitude.
    /// dF/dφ = 1/√(1 - k²sin²φ)
    pub fn deriv_amplitude(&self) -> String {
        format!(
            "1/√(1 - {}² * sin²({}))",
            self.modulus, self.amplitude
        )
    }

    /// Returns the derivative with respect to modulus.
    pub fn deriv_modulus(&self) -> String {
        format!(
            "E({0}, {1})/(({1})*(1-{1}²)) - F({0}, {1})/{1}",
            self.amplitude, self.modulus
        )
    }
}

impl std::fmt::Display for EllipticF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_complete {
            write!(f, "K({})", self.modulus)
        } else {
            write!(f, "F({}, {})", self.amplitude, self.modulus)
        }
    }
}

/// Incomplete elliptic integral of the second kind E(φ, k).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EllipticE {
    /// Amplitude φ.
    pub amplitude: String,
    /// Modulus k.
    pub modulus: String,
    /// Whether this is complete.
    pub is_complete: bool,
}

impl EllipticE {
    /// Creates E(φ, k).
    pub fn incomplete(amplitude: &str, modulus: &str) -> Self {
        Self {
            amplitude: amplitude.to_string(),
            modulus: modulus.to_string(),
            is_complete: false,
        }
    }

    /// Creates the complete elliptic integral E(k).
    pub fn complete(modulus: &str) -> Self {
        Self {
            amplitude: "π/2".to_string(),
            modulus: modulus.to_string(),
            is_complete: true,
        }
    }

    /// Returns the derivative with respect to amplitude.
    /// dE/dφ = √(1 - k²sin²φ)
    pub fn deriv_amplitude(&self) -> String {
        format!(
            "√(1 - {}² * sin²({}))",
            self.modulus, self.amplitude
        )
    }
}

impl std::fmt::Display for EllipticE {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_complete {
            write!(f, "E({})", self.modulus)
        } else {
            write!(f, "E({}, {})", self.amplitude, self.modulus)
        }
    }
}

/// Incomplete elliptic integral of the third kind Π(n; φ, k).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EllipticPi {
    /// Characteristic n.
    pub characteristic: String,
    /// Amplitude φ.
    pub amplitude: String,
    /// Modulus k.
    pub modulus: String,
    /// Whether this is complete.
    pub is_complete: bool,
}

impl EllipticPi {
    /// Creates Π(n; φ, k).
    pub fn incomplete(characteristic: &str, amplitude: &str, modulus: &str) -> Self {
        Self {
            characteristic: characteristic.to_string(),
            amplitude: amplitude.to_string(),
            modulus: modulus.to_string(),
            is_complete: false,
        }
    }

    /// Creates the complete elliptic integral Π(n; k).
    pub fn complete(characteristic: &str, modulus: &str) -> Self {
        Self {
            characteristic: characteristic.to_string(),
            amplitude: "π/2".to_string(),
            modulus: modulus.to_string(),
            is_complete: true,
        }
    }
}

impl std::fmt::Display for EllipticPi {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_complete {
            write!(f, "Π({}; {})", self.characteristic, self.modulus)
        } else {
            write!(
                f, "Π({}; {}, {})",
                self.characteristic, self.amplitude, self.modulus
            )
        }
    }
}

/// Known special values of complete elliptic integrals.
pub struct EllipticValues;

impl EllipticValues {
    /// K(0) = π/2
    pub fn k_at_zero() -> String {
        "π/2".to_string()
    }

    /// E(0) = π/2
    pub fn e_at_zero() -> String {
        "π/2".to_string()
    }

    /// K(1) = ∞ (singular)
    pub fn k_at_one() -> String {
        "∞".to_string()
    }

    /// E(1) = 1
    pub fn e_at_one() -> String {
        "1".to_string()
    }

    /// Legendre's relation: E(k)K(k') + E(k')K(k) - K(k)K(k') = π/2
    /// where k' = √(1-k²)
    pub fn legendre_relation() -> String {
        "E(k)·K(k') + E(k')·K(k) - K(k)·K(k') = π/2, where k' = √(1-k²)".to_string()
    }
}

/// Integrals that reduce to elliptic functions.
pub struct EllipticReductions;

impl EllipticReductions {
    /// ∫ 1/√(1-x⁴) dx = F(arcsin(x), i) / √2
    pub fn quartic_sqrt_reciprocal() -> String {
        "F(arcsin(x), i) / √2".to_string()
    }

    /// ∫ √(1-k²x²)/√(1-x²) dx = E(arcsin(x), k)
    pub fn standard_e_form() -> String {
        "E(arcsin(x), k)".to_string()
    }

    /// Arc length of ellipse with semi-axes a, b
    pub fn ellipse_arc_length() -> String {
        "4a · E(k) where k = √(1 - b²/a²)".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_elliptic_f_incomplete() {
        let f = EllipticF::incomplete("φ", "k");
        assert_eq!(format!("{}", f), "F(φ, k)");
        assert!(!f.is_complete);
    }

    #[test]
    fn test_elliptic_f_complete() {
        let k = EllipticF::complete("k");
        assert_eq!(format!("{}", k), "K(k)");
        assert!(k.is_complete);
    }

    #[test]
    fn test_elliptic_e() {
        let e = EllipticE::incomplete("φ", "k");
        assert_eq!(format!("{}", e), "E(φ, k)");

        let e_complete = EllipticE::complete("k");
        assert_eq!(format!("{}", e_complete), "E(k)");
    }

    #[test]
    fn test_elliptic_pi() {
        let pi = EllipticPi::incomplete("n", "φ", "k");
        assert_eq!(format!("{}", pi), "Π(n; φ, k)");
    }

    #[test]
    fn test_special_values() {
        assert_eq!(EllipticValues::k_at_zero(), "π/2");
        assert_eq!(EllipticValues::e_at_one(), "1");
    }
}
