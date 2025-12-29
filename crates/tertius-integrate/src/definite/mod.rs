//! Definite Integral Evaluation Framework
//!
//! This module provides functionality for evaluating definite integrals
//! using antiderivatives with proper limit handling for singularities
//! and infinite bounds.
//!
//! # Example
//!
//! ```ignore
//! use tertius_integrate::definite::{definite_integral, Bound};
//!
//! // Evaluate ∫₀^∞ e^(-x) dx = 1
//! let result = definite_integral_numerical(
//!     &|x| (-x).exp(),
//!     Bound::Finite(0.0),
//!     Bound::PosInfinity,
//!     1e-10,
//!     1000,
//! );
//! ```

pub mod evaluate;
pub mod improper;

use crate::numerical::{
    adaptive_integrate, adaptive_integrate_to_infinity, adaptive_integrate_from_neg_infinity,
    adaptive_integrate_full_line, AdaptiveResult,
};

/// Represents an integration bound.
#[derive(Clone, Debug, PartialEq)]
pub enum Bound {
    /// A finite bound value
    Finite(f64),
    /// Positive infinity (+∞)
    PosInfinity,
    /// Negative infinity (-∞)
    NegInfinity,
}

/// Result of a definite integral computation.
#[derive(Clone, Debug)]
pub struct DefiniteResult {
    /// Computed integral value
    pub value: f64,
    /// Estimated error
    pub error: f64,
    /// Number of function evaluations
    pub evaluations: usize,
    /// Whether the computation converged
    pub converged: bool,
    /// Whether the integral is proper (no singularities)
    pub is_proper: bool,
}

/// Evaluates a definite integral numerically.
///
/// This function automatically handles:
/// - Finite intervals
/// - Semi-infinite intervals [a, ∞) or (-∞, b]
/// - Full line integrals (-∞, ∞)
/// - Interior singularities (via adaptive subdivision)
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `lower` - Lower bound
/// * `upper` - Upper bound
/// * `tolerance` - Desired accuracy (absolute tolerance)
/// * `max_subdivisions` - Maximum number of adaptive subdivisions
///
/// # Example
///
/// ```ignore
/// use tertius_integrate::definite::{definite_integral_numerical, Bound};
///
/// // ∫₀¹ x² dx = 1/3
/// let result = definite_integral_numerical(
///     &|x| x * x,
///     Bound::Finite(0.0),
///     Bound::Finite(1.0),
///     1e-10,
///     1000,
/// );
/// assert!((result.value - 1.0/3.0).abs() < 1e-9);
/// ```
pub fn definite_integral_numerical<F: Fn(f64) -> f64>(
    f: &F,
    lower: Bound,
    upper: Bound,
    tolerance: f64,
    max_subdivisions: usize,
) -> DefiniteResult {
    match (&lower, &upper) {
        (Bound::Finite(a), Bound::Finite(b)) => {
            // Standard finite interval
            let result = adaptive_integrate(f, *a, *b, tolerance, tolerance, max_subdivisions);
            DefiniteResult {
                value: result.value,
                error: result.error,
                evaluations: result.evaluations,
                converged: result.converged,
                is_proper: true,
            }
        }
        (Bound::Finite(a), Bound::PosInfinity) => {
            // Semi-infinite [a, ∞)
            let result = adaptive_integrate_to_infinity(f, *a, tolerance, tolerance, max_subdivisions);
            DefiniteResult {
                value: result.value,
                error: result.error,
                evaluations: result.evaluations,
                converged: result.converged,
                is_proper: false,
            }
        }
        (Bound::NegInfinity, Bound::Finite(b)) => {
            // Semi-infinite (-∞, b]
            let result = adaptive_integrate_from_neg_infinity(f, *b, tolerance, tolerance, max_subdivisions);
            DefiniteResult {
                value: result.value,
                error: result.error,
                evaluations: result.evaluations,
                converged: result.converged,
                is_proper: false,
            }
        }
        (Bound::NegInfinity, Bound::PosInfinity) => {
            // Full line integral
            let result = adaptive_integrate_full_line(f, tolerance, tolerance, max_subdivisions);
            DefiniteResult {
                value: result.value,
                error: result.error,
                evaluations: result.evaluations,
                converged: result.converged,
                is_proper: false,
            }
        }
        (Bound::PosInfinity, Bound::PosInfinity) | (Bound::NegInfinity, Bound::NegInfinity) => {
            // Degenerate case: same infinite bound
            DefiniteResult {
                value: 0.0,
                error: 0.0,
                evaluations: 0,
                converged: true,
                is_proper: false,
            }
        }
        (Bound::PosInfinity, _) | (_, Bound::NegInfinity) => {
            // Invalid bounds (a > b in the infinite sense)
            DefiniteResult {
                value: f64::NAN,
                error: f64::INFINITY,
                evaluations: 0,
                converged: false,
                is_proper: false,
            }
        }
    }
}

/// Evaluates a definite integral with known singularities.
///
/// Splits the integral at singularity points to avoid numerical issues.
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `lower` - Lower bound
/// * `upper` - Upper bound
/// * `singularities` - Points where the integrand has singularities
/// * `tolerance` - Desired accuracy
/// * `max_subdivisions` - Maximum number of adaptive subdivisions
pub fn definite_integral_with_singularities<F: Fn(f64) -> f64>(
    f: &F,
    lower: f64,
    upper: f64,
    singularities: &[f64],
    tolerance: f64,
    max_subdivisions: usize,
) -> DefiniteResult {
    use crate::numerical::adaptive_integrate_with_singularities;

    let result = adaptive_integrate_with_singularities(
        f, lower, upper, singularities, tolerance, tolerance, max_subdivisions
    );

    DefiniteResult {
        value: result.value,
        error: result.error,
        evaluations: result.evaluations,
        converged: result.converged,
        is_proper: singularities.is_empty(),
    }
}

/// Computes the principal value of an integral with a simple pole.
///
/// For integrals of the form ∫_a^b f(x) dx where f has a simple pole at c ∈ (a, b),
/// computes the Cauchy principal value:
///
/// P.V. ∫_a^b f(x) dx = lim_{ε→0} [∫_a^{c-ε} f(x) dx + ∫_{c+ε}^b f(x) dx]
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `pole` - Location of the simple pole
/// * `tolerance` - Desired accuracy
/// * `max_subdivisions` - Maximum number of subdivisions
pub fn principal_value<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    pole: f64,
    tolerance: f64,
    max_subdivisions: usize,
) -> DefiniteResult {
    if pole <= a || pole >= b {
        // Pole is outside interval, just integrate normally
        return definite_integral_numerical(
            f,
            Bound::Finite(a),
            Bound::Finite(b),
            tolerance,
            max_subdivisions,
        );
    }

    // Compute symmetric approach to the pole
    let epsilon = 1e-10;

    // Integrate on each side
    let left = adaptive_integrate(f, a, pole - epsilon, tolerance / 2.0, tolerance / 2.0, max_subdivisions / 2);
    let right = adaptive_integrate(f, pole + epsilon, b, tolerance / 2.0, tolerance / 2.0, max_subdivisions / 2);

    DefiniteResult {
        value: left.value + right.value,
        error: left.error + right.error,
        evaluations: left.evaluations + right.evaluations,
        converged: left.converged && right.converged,
        is_proper: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_finite_integral() {
        // ∫₀¹ x² dx = 1/3
        let result = definite_integral_numerical(
            &|x| x * x,
            Bound::Finite(0.0),
            Bound::Finite(1.0),
            1e-10,
            100,
        );
        assert!((result.value - 1.0 / 3.0).abs() < 1e-10);
        assert!(result.is_proper);
    }

    #[test]
    fn test_semi_infinite_integral() {
        // ∫₀^∞ e^(-x) dx = 1
        let result = definite_integral_numerical(
            &|x| (-x).exp(),
            Bound::Finite(0.0),
            Bound::PosInfinity,
            1e-6,
            200,
        );
        assert!((result.value - 1.0).abs() < 1e-4);
        assert!(!result.is_proper);
    }

    #[test]
    fn test_full_line_integral() {
        // ∫_{-∞}^∞ e^(-x²) dx = √π
        let result = definite_integral_numerical(
            &|x| (-x * x).exp(),
            Bound::NegInfinity,
            Bound::PosInfinity,
            1e-8,
            200,
        );
        assert!((result.value - PI.sqrt()).abs() < 1e-6);
        assert!(!result.is_proper);
    }

    #[test]
    fn test_integral_with_singularity() {
        // ∫₀¹ 1/√x dx = 2 (integrable singularity at 0)
        let result = definite_integral_with_singularities(
            &|x| if x > 1e-12 { 1.0 / x.sqrt() } else { 0.0 },
            1e-10,
            1.0,
            &[0.0],
            1e-4,
            500,
        );
        assert!((result.value - 2.0).abs() < 0.05);
    }

    #[test]
    fn test_principal_value() {
        // P.V. ∫₋₁¹ 1/x dx = 0 (by symmetry)
        let result = principal_value(
            &|x| if x.abs() > 1e-12 { 1.0 / x } else { 0.0 },
            -1.0,
            1.0,
            0.0,
            1e-10,
            100,
        );
        assert!(result.value.abs() < 1e-8);
    }

    #[test]
    fn test_invalid_bounds() {
        // +∞ as lower bound should give NaN
        let result = definite_integral_numerical(
            &|x| x,
            Bound::PosInfinity,
            Bound::Finite(0.0),
            1e-10,
            100,
        );
        assert!(result.value.is_nan());
    }
}
