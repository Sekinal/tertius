//! Special Definite Integral Recognition
//!
//! This module recognizes and evaluates special definite integrals that have
//! known closed-form values, including:
//!
//! - **Dirichlet integrals**: ∫₀^∞ sin(x)/x dx = π/2
//! - **Beta/Gamma integrals**: ∫₀^1 x^(a-1)(1-x)^(b-1) dx = B(a,b)
//! - **Gaussian integrals**: ∫_{-∞}^∞ e^(-x²) dx = √π
//! - **Golden ratio integral**: The special integral involving √((1+x)/(1-x))
//!
//! # Example
//!
//! ```ignore
//! use tertius_integrate::special_definite::{recognize_special_integral, GoldenRatioIntegral};
//!
//! let result = GoldenRatioIntegral::evaluate_exact();
//! println!("Golden ratio integral = {}", result); // 4π·arccot(φ)
//! ```

pub mod golden_ratio;
pub mod dirichlet;
pub mod beta_gamma;

pub use golden_ratio::{GoldenRatioIntegral, GoldenRatioResult};

use std::f64::consts::PI;

/// Special mathematical constants.
pub struct MathConstants;

impl MathConstants {
    /// The golden ratio φ = (1 + √5)/2 ≈ 1.618033988749895
    pub const PHI: f64 = 1.618033988749894848204586834365638117720309179805762862135;

    /// The reciprocal of the golden ratio 1/φ = φ - 1 ≈ 0.618033988749895
    pub const PHI_RECIPROCAL: f64 = 0.6180339887498948482045868343656381177203091798057628621354;

    /// Euler-Mascheroni constant γ ≈ 0.5772156649015329
    pub const EULER_GAMMA: f64 = 0.5772156649015328606065120900824024310421593359399235988057;

    /// Catalan's constant G ≈ 0.9159655941772190
    pub const CATALAN: f64 = 0.9159655941772190150546035149323841107741493742816721342664;
}

/// Result of recognizing a special integral.
#[derive(Clone, Debug)]
pub enum SpecialIntegralResult {
    /// Dirichlet integral pattern
    Dirichlet { value: f64, description: String },
    /// Beta function integral
    Beta { a: f64, b: f64, value: f64 },
    /// Gamma function integral
    Gamma { n: f64, value: f64 },
    /// Gaussian integral
    Gaussian { value: f64 },
    /// Golden ratio integral
    GoldenRatio(GoldenRatioResult),
    /// Numerical result (no special form recognized)
    Numerical { value: f64, error: f64 },
    /// Not recognized
    Unknown,
}

/// Attempts to recognize a special definite integral pattern.
pub fn recognize_special_integral<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
) -> SpecialIntegralResult {
    // Try Gaussian integral pattern
    if let Some(result) = try_gaussian_integral(f, a, b, tolerance) {
        return result;
    }

    // Try Dirichlet-type patterns
    if let Some(result) = try_dirichlet_integral(f, a, b, tolerance) {
        return result;
    }

    // Fallback to numerical
    use crate::numerical::adaptive_integrate;
    let result = adaptive_integrate(f, a, b, tolerance, tolerance, 1000);
    SpecialIntegralResult::Numerical {
        value: result.value,
        error: result.error,
    }
}

/// Checks if the integrand matches a Gaussian pattern.
fn try_gaussian_integral<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
) -> Option<SpecialIntegralResult> {
    // Check if integrand looks like e^(-x²)
    // Sample at a few points
    let samples = [0.0, 0.5, 1.0, 1.5, 2.0];
    let expected: Vec<f64> = samples.iter().map(|&x: &f64| (-x * x).exp()).collect();
    let actual: Vec<f64> = samples.iter().map(|&x| f(x)).collect();

    let is_gaussian = expected
        .iter()
        .zip(actual.iter())
        .all(|(e, a)| (e - a).abs() < tolerance * 10.0);

    if is_gaussian {
        // Check bounds
        if a.is_infinite() && a < 0.0 && b.is_infinite() && b > 0.0 {
            return Some(SpecialIntegralResult::Gaussian { value: PI.sqrt() });
        }
    }

    None
}

/// Checks if the integrand matches a Dirichlet pattern.
fn try_dirichlet_integral<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
) -> Option<SpecialIntegralResult> {
    // Check if integrand looks like sin(x)/x
    if (a - 0.0).abs() < tolerance && b.is_infinite() && b > 0.0 {
        // Check at small x (should approach 1)
        let at_small = f(0.001);
        if (at_small - 1.0).abs() < 0.1 {
            // Verify oscillation pattern
            let at_pi = f(PI);
            let at_2pi = f(2.0 * PI);
            if at_pi.abs() < 0.1 && at_2pi.abs() < 0.1 {
                return Some(SpecialIntegralResult::Dirichlet {
                    value: PI / 2.0,
                    description: "∫₀^∞ sin(x)/x dx = π/2".to_string(),
                });
            }
        }
    }

    None
}

/// Computes the Beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b).
pub fn beta(a: f64, b: f64) -> f64 {
    use std::f64::consts::PI;

    // Use the relation B(a,b) = ∫₀^1 t^(a-1)(1-t)^(b-1) dt
    // For integer/half-integer arguments, use exact formulas

    // General case via gamma functions
    gamma(a) * gamma(b) / gamma(a + b)
}

/// Computes the Gamma function Γ(x).
///
/// Uses Lanczos approximation for general values.
pub fn gamma(x: f64) -> f64 {
    if x <= 0.0 && x == x.floor() {
        return f64::INFINITY; // Poles at non-positive integers
    }

    if x < 0.5 {
        // Reflection formula: Γ(1-x)Γ(x) = π/sin(πx)
        PI / ((PI * x).sin() * gamma(1.0 - x))
    } else {
        // Lanczos approximation
        lanczos_gamma(x)
    }
}

/// Lanczos approximation for Γ(x) with x ≥ 0.5.
fn lanczos_gamma(x: f64) -> f64 {
    // Lanczos coefficients for g=7
    const G: f64 = 7.0;
    const COEFFS: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    let x = x - 1.0;
    let mut ag = COEFFS[0];
    for i in 1..9 {
        ag += COEFFS[i] / (x + i as f64);
    }

    let t = x + G + 0.5;
    (2.0 * PI).sqrt() * t.powf(x + 0.5) * (-t).exp() * ag
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phi_constants() {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((MathConstants::PHI - phi).abs() < 1e-15);
        assert!((MathConstants::PHI_RECIPROCAL - (phi - 1.0)).abs() < 1e-15);
    }

    #[test]
    fn test_gamma_integer() {
        // Γ(n) = (n-1)! for positive integers
        assert!((gamma(1.0) - 1.0).abs() < 1e-10); // 0!
        assert!((gamma(2.0) - 1.0).abs() < 1e-10); // 1!
        assert!((gamma(3.0) - 2.0).abs() < 1e-10); // 2!
        assert!((gamma(4.0) - 6.0).abs() < 1e-10); // 3!
        assert!((gamma(5.0) - 24.0).abs() < 1e-10); // 4!
    }

    #[test]
    fn test_gamma_half() {
        // Γ(1/2) = √π
        assert!((gamma(0.5) - PI.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_beta_symmetry() {
        // B(a, b) = B(b, a)
        assert!((beta(2.0, 3.0) - beta(3.0, 2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_beta_integer() {
        // B(a, b) = (a-1)!(b-1)!/(a+b-1)! for integers
        // B(2, 3) = 1!2!/4! = 2/24 = 1/12
        assert!((beta(2.0, 3.0) - 1.0 / 12.0).abs() < 1e-10);
    }
}
