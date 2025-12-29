//! Algebraic Function Integration
//!
//! This module implements integration of algebraic functions, i.e.,
//! functions involving radicals like √P(x) where P is a polynomial.
//!
//! # Theory
//!
//! An algebraic function is a function y = y(x) satisfying a polynomial
//! equation P(x, y) = 0. For square roots, y = √Q(x) satisfies y² - Q(x) = 0.
//!
//! The **genus** of the algebraic curve determines the integration method:
//! - Genus 0: Rationalizable (e.g., √(ax² + bx + c)) → elementary
//! - Genus 1: Elliptic integrals (cubic/quartic radicands)
//! - Genus ≥ 2: Hyperelliptic/Abelian integrals
//!
//! # Implemented Algorithms
//!
//! - **Rationalization** for genus 0 cases
//! - **Davenport-Trager** for general algebraic functions
//! - **Elliptic reduction** to Legendre normal form
//!
//! # Example
//!
//! ```ignore
//! use tertius_integrate::algebraic::{AlgebraicIntegrand, integrate_algebraic};
//!
//! // Integrate √(1-x²)
//! let integrand = AlgebraicIntegrand::sqrt_of_poly(vec![-1, 0, 1]); // 1 - x²
//! let result = integrate_algebraic(&integrand);
//! ```

pub mod rationalize;
pub mod sqrt_rational;
pub mod elliptic;

use std::fmt;

/// Represents the genus of an algebraic curve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Genus {
    /// Genus 0: Can be rationalized
    Rational,
    /// Genus 1: Elliptic curve (cubic or quartic under square root)
    Elliptic,
    /// Genus ≥ 2: Hyperelliptic or general algebraic curve
    Hyperelliptic(usize),
}

impl Genus {
    /// Computes the genus from a polynomial degree under a square root.
    ///
    /// For y² = P(x) with deg(P) = n:
    /// - n ≤ 2: genus 0
    /// - n = 3, 4: genus 1
    /// - n ≥ 5: genus = floor((n-1)/2)
    pub fn from_degree(n: usize) -> Self {
        match n {
            0 | 1 | 2 => Genus::Rational,
            3 | 4 => Genus::Elliptic,
            n => Genus::Hyperelliptic((n - 1) / 2),
        }
    }
}

/// An integrand involving a square root of a rational function.
///
/// Represents f(x) = R(x, √(P(x)/Q(x))) where R is rational in both arguments.
#[derive(Clone, Debug)]
pub struct AlgebraicIntegrand {
    /// Polynomial under the square root (numerator)
    pub radicand_num: Vec<f64>,
    /// Polynomial under the square root (denominator), or [1.0] for polynomial
    pub radicand_den: Vec<f64>,
    /// Coefficients of the rational function outside the radical
    /// Represented as (numerator_poly, denominator_poly)
    pub rational_part: (Vec<f64>, Vec<f64>),
    /// Whether the square root appears in the numerator (true) or denominator (false)
    pub sqrt_in_numerator: bool,
}

impl AlgebraicIntegrand {
    /// Creates an integrand √P(x) (just the square root).
    pub fn sqrt_of_poly(coeffs: Vec<f64>) -> Self {
        Self {
            radicand_num: coeffs,
            radicand_den: vec![1.0],
            rational_part: (vec![1.0], vec![1.0]),
            sqrt_in_numerator: true,
        }
    }

    /// Creates an integrand R(x) * √P(x).
    pub fn rational_times_sqrt(r_num: Vec<f64>, r_den: Vec<f64>, radicand: Vec<f64>) -> Self {
        Self {
            radicand_num: radicand,
            radicand_den: vec![1.0],
            rational_part: (r_num, r_den),
            sqrt_in_numerator: true,
        }
    }

    /// Creates an integrand R(x) / √P(x).
    pub fn rational_over_sqrt(r_num: Vec<f64>, r_den: Vec<f64>, radicand: Vec<f64>) -> Self {
        Self {
            radicand_num: radicand,
            radicand_den: vec![1.0],
            rational_part: (r_num, r_den),
            sqrt_in_numerator: false,
        }
    }

    /// Creates the integrand for √((1+x)/(1-x)).
    pub fn sqrt_ratio() -> Self {
        Self {
            radicand_num: vec![1.0, 1.0],    // 1 + x
            radicand_den: vec![1.0, -1.0],   // 1 - x
            rational_part: (vec![1.0], vec![1.0]),
            sqrt_in_numerator: true,
        }
    }

    /// Returns the degree of the radicand polynomial.
    pub fn radicand_degree(&self) -> usize {
        degree(&self.radicand_num).max(degree(&self.radicand_den))
    }

    /// Determines the genus of the associated algebraic curve.
    pub fn genus(&self) -> Genus {
        // For √(P/Q), we look at the degree of P*Q (after clearing squares)
        let combined_deg = degree(&self.radicand_num) + degree(&self.radicand_den);
        Genus::from_degree(combined_deg)
    }

    /// Evaluates the integrand at a point.
    pub fn evaluate(&self, x: f64) -> f64 {
        let p_val = eval_poly(&self.radicand_num, x);
        let q_val = eval_poly(&self.radicand_den, x);
        let sqrt_val = (p_val / q_val).sqrt();

        let r_num = eval_poly(&self.rational_part.0, x);
        let r_den = eval_poly(&self.rational_part.1, x);

        if self.sqrt_in_numerator {
            (r_num / r_den) * sqrt_val
        } else {
            (r_num / r_den) / sqrt_val
        }
    }
}

/// Result of integrating an algebraic function.
#[derive(Clone, Debug)]
pub enum AlgebraicIntegralResult {
    /// Elementary antiderivative found
    Elementary {
        /// Rational part of the antiderivative
        rational: String,
        /// Logarithmic part (sum of c_i * ln(arg_i))
        logs: Vec<(f64, String)>,
        /// Arctangent part
        arctans: Vec<(f64, String)>,
    },
    /// Elliptic integral (first, second, or third kind)
    Elliptic {
        /// Kind of elliptic integral (1, 2, or 3)
        kind: u8,
        /// Amplitude φ
        amplitude: String,
        /// Modulus k
        modulus: f64,
        /// Parameter n (for third kind only)
        parameter: Option<f64>,
    },
    /// Hyperelliptic or Abelian integral (no elementary form)
    NonElementary {
        /// Description of the integral type
        description: String,
        /// Genus of the associated curve
        genus: usize,
    },
    /// Numerical value only (for definite integrals)
    Numerical(f64),
}

impl fmt::Display for AlgebraicIntegralResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlgebraicIntegralResult::Elementary { rational, logs, arctans } => {
                write!(f, "{}", rational)?;
                for (coef, arg) in logs {
                    if *coef >= 0.0 {
                        write!(f, " + {}*ln({})", coef, arg)?;
                    } else {
                        write!(f, " - {}*ln({})", -coef, arg)?;
                    }
                }
                for (coef, arg) in arctans {
                    if *coef >= 0.0 {
                        write!(f, " + {}*arctan({})", coef, arg)?;
                    } else {
                        write!(f, " - {}*arctan({})", -coef, arg)?;
                    }
                }
                Ok(())
            }
            AlgebraicIntegralResult::Elliptic { kind, amplitude, modulus, parameter } => {
                match kind {
                    1 => write!(f, "F({}, {})", amplitude, modulus),
                    2 => write!(f, "E({}, {})", amplitude, modulus),
                    3 => write!(f, "Π({}; {}, {})", parameter.unwrap_or(0.0), amplitude, modulus),
                    _ => write!(f, "Elliptic[{}]", kind),
                }
            }
            AlgebraicIntegralResult::NonElementary { description, genus } => {
                write!(f, "{} (genus {})", description, genus)
            }
            AlgebraicIntegralResult::Numerical(val) => {
                write!(f, "{}", val)
            }
        }
    }
}

/// Integrates an algebraic function.
///
/// Attempts to find an elementary antiderivative; if not possible,
/// expresses the result in terms of elliptic or hyperelliptic integrals.
pub fn integrate_algebraic(integrand: &AlgebraicIntegrand) -> AlgebraicIntegralResult {
    match integrand.genus() {
        Genus::Rational => {
            // Try to rationalize and integrate
            rationalize::integrate_rationalizable(integrand)
        }
        Genus::Elliptic => {
            // Reduce to elliptic integrals
            elliptic::integrate_elliptic(integrand)
        }
        Genus::Hyperelliptic(g) => {
            // Cannot express elementarily
            AlgebraicIntegralResult::NonElementary {
                description: format!("Hyperelliptic integral"),
                genus: g,
            }
        }
    }
}

/// Numerically integrates an algebraic function over [a, b].
pub fn integrate_algebraic_numerical(
    integrand: &AlgebraicIntegrand,
    a: f64,
    b: f64,
    tolerance: f64,
) -> f64 {
    use crate::numerical::adaptive_integrate;

    let f = |x: f64| integrand.evaluate(x);
    let result = adaptive_integrate(&f, a, b, tolerance, tolerance, 1000);
    result.value
}

// Helper functions for polynomial operations

fn degree(p: &[f64]) -> usize {
    if p.is_empty() {
        return 0;
    }
    for i in (0..p.len()).rev() {
        if p[i].abs() > 1e-14 {
            return i;
        }
    }
    0
}

fn eval_poly(p: &[f64], x: f64) -> f64 {
    let mut result = 0.0;
    let mut power = 1.0;
    for &coef in p {
        result += coef * power;
        power *= x;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genus_classification() {
        // √(ax + b) → genus 0
        assert_eq!(Genus::from_degree(1), Genus::Rational);

        // √(ax² + bx + c) → genus 0
        assert_eq!(Genus::from_degree(2), Genus::Rational);

        // √(x³ + ...) → genus 1 (elliptic)
        assert_eq!(Genus::from_degree(3), Genus::Elliptic);

        // √(x⁴ + ...) → genus 1 (elliptic)
        assert_eq!(Genus::from_degree(4), Genus::Elliptic);

        // √(x⁵ + ...) → genus 2 (hyperelliptic)
        assert_eq!(Genus::from_degree(5), Genus::Hyperelliptic(2));
    }

    #[test]
    fn test_integrand_evaluation() {
        // √(1 - x²) at x = 0 should be 1
        let integrand = AlgebraicIntegrand::sqrt_of_poly(vec![1.0, 0.0, -1.0]);
        assert!((integrand.evaluate(0.0) - 1.0).abs() < 1e-10);

        // √(1 - x²) at x = 0.5 should be √0.75
        let expected = 0.75_f64.sqrt();
        assert!((integrand.evaluate(0.5) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sqrt_ratio() {
        // √((1+x)/(1-x)) at x = 0 should be 1
        let integrand = AlgebraicIntegrand::sqrt_ratio();
        assert!((integrand.evaluate(0.0) - 1.0).abs() < 1e-10);

        // √((1+0.5)/(1-0.5)) = √3 at x = 0.5
        let expected = 3.0_f64.sqrt();
        assert!((integrand.evaluate(0.5) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_numerical_integration() {
        // ∫₀^1 √(1-x²) dx = π/4
        let integrand = AlgebraicIntegrand::sqrt_of_poly(vec![1.0, 0.0, -1.0]);
        let result = integrate_algebraic_numerical(&integrand, 0.0, 1.0, 1e-8);
        let expected = std::f64::consts::PI / 4.0;
        assert!((result - expected).abs() < 1e-6);
    }
}
