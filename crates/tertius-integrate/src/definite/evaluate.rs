//! Antiderivative Evaluation at Bounds
//!
//! This module provides utilities for evaluating antiderivatives at
//! integration bounds, handling limits at singularities and infinity.

use crate::numerical::AdaptiveResult;

/// Result of evaluating an antiderivative at bounds.
#[derive(Clone, Debug)]
pub struct BoundEvaluationResult {
    /// Value at upper bound: F(b)
    pub upper_value: Option<f64>,
    /// Value at lower bound: F(a)
    pub lower_value: Option<f64>,
    /// Final integral value: F(b) - F(a)
    pub integral_value: f64,
    /// Whether limits were needed (discontinuity or infinity)
    pub needed_limits: bool,
    /// Description of any issues encountered
    pub notes: Vec<String>,
}

/// Evaluates a function at a point, with handling for singularities.
///
/// If direct evaluation fails or returns NaN/Infinity, attempts to
/// compute a limit from one or both sides.
///
/// # Arguments
///
/// * `f` - The function to evaluate
/// * `x` - The point at which to evaluate
/// * `from_left` - Whether to approach from the left (x - ε) or right (x + ε)
pub fn evaluate_with_limit<F: Fn(f64) -> f64>(
    f: &F,
    x: f64,
    from_left: bool,
) -> Option<f64> {
    // First try direct evaluation
    let direct = f(x);
    if direct.is_finite() {
        return Some(direct);
    }

    // Try approaching from the specified direction
    let epsilon_start = 1e-2;
    let epsilon_min = 1e-14;

    let mut epsilon = epsilon_start;
    let mut prev_value = f64::NAN;
    let mut converging = false;

    while epsilon > epsilon_min {
        let point = if from_left { x - epsilon } else { x + epsilon };
        let value = f(point);

        if value.is_finite() {
            if prev_value.is_finite() {
                // Check if we're converging
                let diff = (value - prev_value).abs();
                if diff < epsilon.sqrt() {
                    converging = true;
                    if diff < 1e-10 {
                        return Some(value);
                    }
                }
            }
            prev_value = value;
        }

        epsilon /= 10.0;
    }

    if converging && prev_value.is_finite() {
        Some(prev_value)
    } else {
        None
    }
}

/// Evaluates an antiderivative F at both bounds to compute ∫_a^b.
///
/// Handles the case where F(a) or F(b) requires a limit computation.
///
/// # Arguments
///
/// * `antiderivative` - The antiderivative F
/// * `a` - Lower bound
/// * `b` - Upper bound
pub fn evaluate_antiderivative_at_bounds<F: Fn(f64) -> f64>(
    antiderivative: &F,
    a: f64,
    b: f64,
) -> BoundEvaluationResult {
    let mut notes = Vec::new();
    let mut needed_limits = false;

    // Evaluate at upper bound
    let upper_direct = antiderivative(b);
    let upper_value = if upper_direct.is_finite() {
        Some(upper_direct)
    } else {
        needed_limits = true;
        notes.push(format!("F(b={}) required limit evaluation", b));
        // Approach from left since we're at upper bound
        evaluate_with_limit(antiderivative, b, true)
    };

    // Evaluate at lower bound
    let lower_direct = antiderivative(a);
    let lower_value = if lower_direct.is_finite() {
        Some(lower_direct)
    } else {
        needed_limits = true;
        notes.push(format!("F(a={}) required limit evaluation", a));
        // Approach from right since we're at lower bound
        evaluate_with_limit(antiderivative, a, false)
    };

    // Compute the integral value
    let integral_value = match (upper_value, lower_value) {
        (Some(fb), Some(fa)) => fb - fa,
        (None, Some(_)) => {
            notes.push("Upper bound limit diverges or DNE".to_string());
            f64::NAN
        }
        (Some(_), None) => {
            notes.push("Lower bound limit diverges or DNE".to_string());
            f64::NAN
        }
        (None, None) => {
            notes.push("Both bounds require limits that diverge or DNE".to_string());
            f64::NAN
        }
    };

    BoundEvaluationResult {
        upper_value,
        lower_value,
        integral_value,
        needed_limits,
        notes,
    }
}

/// Evaluates an antiderivative at infinity.
///
/// Computes lim_{x→∞} F(x) or lim_{x→-∞} F(x).
///
/// # Arguments
///
/// * `antiderivative` - The antiderivative F
/// * `positive_infinity` - Whether to compute lim as x→+∞ (true) or x→-∞ (false)
pub fn evaluate_at_infinity<F: Fn(f64) -> f64>(
    antiderivative: &F,
    positive_infinity: bool,
) -> Option<f64> {
    // Try increasingly large values and check for convergence
    let start = if positive_infinity { 1.0 } else { -1.0 };
    let mut x = start;
    let mut prev_value = antiderivative(x);

    if !prev_value.is_finite() {
        return None;
    }

    for _ in 0..20 {
        x *= 10.0; // Go to ±10, ±100, ±1000, etc.
        let value = antiderivative(x);

        if !value.is_finite() {
            // Diverging
            return None;
        }

        // Check convergence
        let diff = (value - prev_value).abs();
        if diff < 1e-10 {
            return Some(value);
        }

        prev_value = value;
    }

    None // Didn't converge
}

/// Computes ∫_a^∞ using antiderivative F.
///
/// Returns F(∞) - F(a).
pub fn evaluate_to_infinity<F: Fn(f64) -> f64>(
    antiderivative: &F,
    a: f64,
) -> BoundEvaluationResult {
    let mut notes = Vec::new();
    let needed_limits = true;

    // Evaluate at lower bound
    let lower_direct = antiderivative(a);
    let lower_value = if lower_direct.is_finite() {
        Some(lower_direct)
    } else {
        notes.push(format!("F(a={}) required limit evaluation", a));
        evaluate_with_limit(antiderivative, a, false)
    };

    // Evaluate limit at infinity
    let upper_value = evaluate_at_infinity(antiderivative, true);
    if upper_value.is_none() {
        notes.push("F(∞) limit diverges or DNE".to_string());
    }

    let integral_value = match (upper_value, lower_value) {
        (Some(finf), Some(fa)) => finf - fa,
        _ => f64::NAN,
    };

    BoundEvaluationResult {
        upper_value,
        lower_value,
        integral_value,
        needed_limits,
        notes,
    }
}

/// Computes ∫_{-∞}^b using antiderivative F.
///
/// Returns F(b) - F(-∞).
pub fn evaluate_from_neg_infinity<F: Fn(f64) -> f64>(
    antiderivative: &F,
    b: f64,
) -> BoundEvaluationResult {
    let mut notes = Vec::new();
    let needed_limits = true;

    // Evaluate at upper bound
    let upper_direct = antiderivative(b);
    let upper_value = if upper_direct.is_finite() {
        Some(upper_direct)
    } else {
        notes.push(format!("F(b={}) required limit evaluation", b));
        evaluate_with_limit(antiderivative, b, true)
    };

    // Evaluate limit at negative infinity
    let lower_value = evaluate_at_infinity(antiderivative, false);
    if lower_value.is_none() {
        notes.push("F(-∞) limit diverges or DNE".to_string());
    }

    let integral_value = match (upper_value, lower_value) {
        (Some(fb), Some(fneginf)) => fb - fneginf,
        _ => f64::NAN,
    };

    BoundEvaluationResult {
        upper_value,
        lower_value,
        integral_value,
        needed_limits,
        notes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_simple_evaluation() {
        // F(x) = x³/3, evaluate at [0, 1]
        let f = |x: f64| x.powi(3) / 3.0;
        let result = evaluate_antiderivative_at_bounds(&f, 0.0, 1.0);

        assert!((result.integral_value - 1.0 / 3.0).abs() < 1e-14);
        assert!(!result.needed_limits);
    }

    #[test]
    fn test_evaluation_with_singularity() {
        // F(x) = ln(x), evaluate at [1, e]
        // ln(1) = 0, ln(e) = 1
        let f = |x: f64| x.ln();
        let result = evaluate_antiderivative_at_bounds(&f, 1.0, std::f64::consts::E);

        assert!((result.integral_value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_evaluation_limit_at_zero() {
        // F(x) = x * ln(x) - x, antiderivative of ln(x)
        // lim_{x→0+} x*ln(x) - x = 0
        let f = |x: f64| {
            if x > 0.0 {
                x * x.ln() - x
            } else {
                f64::NEG_INFINITY
            }
        };
        let result = evaluate_antiderivative_at_bounds(&f, 0.0, 1.0);

        // F(1) - F(0+) = (0 - 1) - 0 = -1
        // This corresponds to ∫₀¹ ln(x) dx = -1
        assert!((result.integral_value - (-1.0)).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_to_infinity() {
        // F(x) = -e^(-x), antiderivative of e^(-x)
        // F(∞) = 0, F(0) = -1, so ∫₀^∞ e^(-x) dx = 0 - (-1) = 1
        let f = |x: f64| -(-x).exp();
        let result = evaluate_to_infinity(&f, 0.0);

        assert!((result.integral_value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_from_neg_infinity() {
        // F(x) = e^x, antiderivative of e^x
        // F(0) = 1, F(-∞) = 0, so ∫_{-∞}^0 e^x dx = 1 - 0 = 1
        let f = |x: f64| x.exp();
        let result = evaluate_from_neg_infinity(&f, 0.0);

        assert!((result.integral_value - 1.0).abs() < 1e-6);
    }
}
