//! Improper Integral Handling
//!
//! This module provides specialized handling for improper integrals:
//! - Integrals with infinite bounds
//! - Integrals with singularities at endpoints or interior points
//! - Oscillatory integrals
//!
//! # Types of Improper Integrals
//!
//! ## Type I: Infinite Bounds
//! - ∫_a^∞ f(x) dx
//! - ∫_{-∞}^b f(x) dx
//! - ∫_{-∞}^∞ f(x) dx
//!
//! ## Type II: Singularities
//! - f has a vertical asymptote at endpoint
//! - f has a vertical asymptote in the interior
//!
//! ## Type III: Oscillatory
//! - ∫_0^∞ sin(x)/x dx (Dirichlet integral)

use crate::numerical::{
    adaptive_integrate, AdaptiveResult,
};

/// Classification of singularity behavior.
#[derive(Clone, Debug, PartialEq)]
pub enum SingularityType {
    /// Integrable singularity (e.g., 1/√x at x=0)
    Integrable,
    /// Simple pole (needs principal value)
    SimplePole,
    /// Higher-order pole (diverges)
    HigherOrderPole,
    /// Oscillatory behavior
    Oscillatory,
    /// Unknown/other
    Unknown,
}

/// Information about a singularity.
#[derive(Clone, Debug)]
pub struct SingularityInfo {
    /// Location of the singularity
    pub location: f64,
    /// Type of singularity
    pub singularity_type: SingularityType,
    /// Estimated order (for poles)
    pub order: Option<f64>,
}

/// Classifies the singularity of f at point x.
///
/// Uses numerical probing to determine the behavior.
pub fn classify_singularity<F: Fn(f64) -> f64>(
    f: &F,
    x: f64,
    approach_from_left: bool,
) -> SingularityInfo {
    let sign = if approach_from_left { -1.0 } else { 1.0 };

    // Probe at decreasing distances from singularity
    let epsilons = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
    let mut values: Vec<f64> = Vec::new();
    let mut log_values: Vec<f64> = Vec::new();
    let mut log_epsilons: Vec<f64> = Vec::new();

    for &eps in &epsilons {
        let point = x + sign * eps;
        let value = f(point).abs();
        if value.is_finite() && value > 0.0 {
            values.push(value);
            log_values.push(value.ln());
            log_epsilons.push(eps.ln());
        }
    }

    if values.len() < 3 {
        return SingularityInfo {
            location: x,
            singularity_type: SingularityType::Unknown,
            order: None,
        };
    }

    // Check if values are bounded AND not growing significantly (no singularity or removable)
    let max_val = values.iter().cloned().fold(0.0_f64, f64::max);
    let min_val = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let growth_ratio = max_val / min_val.max(1e-15);

    // If all values bounded AND ratio is small, it's effectively integrable
    if max_val < 1e6 && growth_ratio < 100.0 {
        return SingularityInfo {
            location: x,
            singularity_type: SingularityType::Integrable,
            order: Some(0.0),
        };
    }

    // Estimate order from log-log slope: |f(x)| ~ C/|x-a|^α
    // ln|f| ~ ln(C) - α * ln|x-a| = ln(C) - α * ln(ε)
    // So slope of ln|f| vs ln(ε) gives -α
    let n = log_values.len();
    let slope = (log_values[n - 1] - log_values[0]) / (log_epsilons[n - 1] - log_epsilons[0]);
    let alpha = -slope;

    // Use tolerance for floating point comparison
    let singularity_type = if alpha < 0.95 {
        // f ~ 1/x^α with α < 1 is integrable
        SingularityType::Integrable
    } else if alpha >= 0.95 && alpha <= 1.05 {
        // f ~ 1/x is a simple pole (α ≈ 1)
        SingularityType::SimplePole
    } else if alpha > 1.05 {
        // f ~ 1/x^α with α > 1 diverges
        SingularityType::HigherOrderPole
    } else {
        SingularityType::Unknown
    };

    SingularityInfo {
        location: x,
        singularity_type,
        order: Some(alpha),
    }
}

/// Integrates over [a, b] with an integrable singularity at a.
///
/// Uses power-law regularization: f(x) ~ C/|x-a|^α with α < 1.
pub fn integrate_left_singularity<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
    max_subdivisions: usize,
) -> AdaptiveResult {
    // Start integration a small distance from the singularity
    let epsilon = 1e-10 * (b - a);
    adaptive_integrate(f, a + epsilon, b, tolerance, tolerance, max_subdivisions)
}

/// Integrates over [a, b] with an integrable singularity at b.
pub fn integrate_right_singularity<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
    max_subdivisions: usize,
) -> AdaptiveResult {
    let epsilon = 1e-10 * (b - a);
    adaptive_integrate(f, a, b - epsilon, tolerance, tolerance, max_subdivisions)
}

/// Integrates using IMT (Iri-Moriguti-Takahashi) transformation.
///
/// The IMT transformation handles endpoint singularities by using
/// a smooth change of variables that clusters quadrature points
/// near the singularities.
///
/// Transforms ∫_a^b f(x) dx to ∫_0^1 f(φ(t)) φ'(t) dt
/// where φ(t) = a + (b-a) * ψ(t) with ψ designed to have
/// vanishing derivatives at t=0 and t=1.
pub fn integrate_imt<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
    max_subdivisions: usize,
) -> AdaptiveResult {
    // Simple IMT transformation: ψ(t) = 0.5 * (1 - cos(πt))
    // This clusters points near endpoints
    let width = b - a;

    let transformed = |t: f64| {
        if t <= 0.0 || t >= 1.0 {
            return 0.0;
        }
        use std::f64::consts::PI;

        // ψ(t) = 0.5 * (1 - cos(πt))
        let psi = 0.5 * (1.0 - (PI * t).cos());
        // ψ'(t) = 0.5 * π * sin(πt)
        let psi_prime = 0.5 * PI * (PI * t).sin();

        let x = a + width * psi;
        f(x) * width * psi_prime
    };

    // Integrate over [ε, 1-ε] to avoid exact endpoints
    let eps = 1e-12;
    adaptive_integrate(&transformed, eps, 1.0 - eps, tolerance, tolerance, max_subdivisions)
}

/// Integrates using double exponential (tanh-sinh) transformation.
///
/// This is one of the most powerful transformations for handling
/// endpoint singularities, achieving rapid convergence even for
/// severe singularities.
///
/// Uses: x = tanh(π/2 * sinh(t)), dx = (π/2 * cosh(t)) / cosh²(π/2 * sinh(t))
pub fn integrate_tanh_sinh<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
    max_subdivisions: usize,
) -> AdaptiveResult {
    use std::f64::consts::PI;

    let mid = (a + b) / 2.0;
    let half_width = (b - a) / 2.0;

    // Transformed integrand: f(mid + half_width * tanh(π/2 * sinh(t))) * jacobian
    let transformed = |t: f64| {
        let sinh_t = t.sinh();
        let cosh_t = t.cosh();

        let arg = PI / 2.0 * sinh_t;
        let tanh_arg = arg.tanh();
        let sech_arg = 1.0 / arg.cosh();

        let x = mid + half_width * tanh_arg;

        // Check if x is within [a, b]
        if x <= a || x >= b {
            return 0.0;
        }

        let jacobian = half_width * (PI / 2.0) * cosh_t * sech_arg * sech_arg;
        let value = f(x);

        if !value.is_finite() || !jacobian.is_finite() {
            0.0
        } else {
            value * jacobian
        }
    };

    // The transformed integral is over (-∞, ∞) but most contribution is in [-5, 5]
    adaptive_integrate(&transformed, -5.0, 5.0, tolerance, tolerance, max_subdivisions)
}

/// Result of oscillatory integral computation.
#[derive(Clone, Debug)]
pub struct OscillatoryResult {
    /// Computed value
    pub value: f64,
    /// Error estimate
    pub error: f64,
    /// Number of periods integrated
    pub periods: usize,
    /// Whether extrapolation was used
    pub extrapolated: bool,
}

/// Integrates oscillatory integrals using Levin-type methods.
///
/// For integrals of the form ∫_a^b f(x) exp(i*ω*g(x)) dx,
/// uses partial sums of periods combined with convergence acceleration.
pub fn integrate_oscillatory_real<F: Fn(f64) -> f64, G: Fn(f64) -> f64>(
    amplitude: &F,
    oscillator: &G,  // sin or cos of (ω*x)
    a: f64,
    b: f64,
    period: f64,     // 2π/ω
    tolerance: f64,
    max_periods: usize,
) -> OscillatoryResult {
    if b - a < period {
        // Less than one period, just integrate directly
        let result = adaptive_integrate(&|x| amplitude(x) * oscillator(x), a, b, tolerance, tolerance, 100);
        return OscillatoryResult {
            value: result.value,
            error: result.error,
            periods: 0,
            extrapolated: false,
        };
    }

    // Integrate period by period and apply convergence acceleration
    let mut partial_sums: Vec<f64> = Vec::new();
    let mut current_sum = 0.0;
    let mut x = a;
    let mut periods_done = 0;

    while x + period <= b && periods_done < max_periods {
        let result = adaptive_integrate(
            &|t| amplitude(t) * oscillator(t),
            x, x + period,
            tolerance, tolerance, 50
        );
        current_sum += result.value;
        partial_sums.push(current_sum);
        x += period;
        periods_done += 1;
    }

    // Handle remaining partial period
    if x < b {
        let result = adaptive_integrate(
            &|t| amplitude(t) * oscillator(t),
            x, b,
            tolerance, tolerance, 50
        );
        current_sum += result.value;
        partial_sums.push(current_sum);
    }

    // Apply Wynn's epsilon algorithm for convergence acceleration
    if partial_sums.len() >= 3 {
        let accelerated = wynn_epsilon(&partial_sums);
        OscillatoryResult {
            value: accelerated,
            error: (accelerated - current_sum).abs(),
            periods: periods_done,
            extrapolated: true,
        }
    } else {
        OscillatoryResult {
            value: current_sum,
            error: tolerance * periods_done as f64,
            periods: periods_done,
            extrapolated: false,
        }
    }
}

/// Wynn's epsilon algorithm for convergence acceleration.
///
/// Given partial sums S_0, S_1, S_2, ..., computes a better
/// estimate of the limit using the epsilon table.
fn wynn_epsilon(partial_sums: &[f64]) -> f64 {
    let n = partial_sums.len();
    if n < 3 {
        return *partial_sums.last().unwrap_or(&0.0);
    }

    // Build epsilon table
    // ε_{-1}(S_n) = 0
    // ε_0(S_n) = S_n
    // ε_{k+1}(S_n) = ε_{k-1}(S_{n+1}) + 1/(ε_k(S_{n+1}) - ε_k(S_n))

    let mut eps_prev: Vec<f64> = vec![0.0; n]; // ε_{k-1}
    let mut eps_curr: Vec<f64> = partial_sums.to_vec(); // ε_k

    let max_iters = (n - 1) / 2;

    for _ in 0..max_iters {
        let mut eps_next = vec![0.0; eps_curr.len() - 1];

        for i in 0..eps_next.len() {
            let diff = eps_curr[i + 1] - eps_curr[i];
            if diff.abs() < 1e-15 {
                eps_next[i] = eps_prev[i + 1];
            } else {
                eps_next[i] = eps_prev[i + 1] + 1.0 / diff;
            }
        }

        if eps_next.is_empty() {
            break;
        }

        eps_prev = eps_curr;
        eps_curr = eps_next;
    }

    // Return the best approximation (even-indexed epsilon values converge to limit)
    *eps_curr.last().unwrap_or(partial_sums.last().unwrap_or(&0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_classify_integrable_singularity() {
        // f(x) = 1/√x has integrable singularity at 0
        let f = |x: f64| 1.0 / x.sqrt();
        let info = classify_singularity(&f, 0.0, false);

        assert_eq!(info.singularity_type, SingularityType::Integrable);
        assert!(info.order.unwrap() < 1.0);
    }

    #[test]
    fn test_classify_simple_pole() {
        // f(x) = 1/x has simple pole at 0
        let f = |x: f64| 1.0 / x;
        let info = classify_singularity(&f, 0.0, false);

        assert_eq!(info.singularity_type, SingularityType::SimplePole);
        assert!((info.order.unwrap() - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_imt_transformation() {
        // ∫₀¹ 1/√x dx = 2
        let f = |x: f64| if x > 0.0 { 1.0 / x.sqrt() } else { 0.0 };
        let result = integrate_imt(&f, 0.0, 1.0, 1e-4, 200);

        assert!((result.value - 2.0).abs() < 0.1);
    }

    #[test]
    fn test_tanh_sinh_transformation() {
        // ∫₀¹ 1/√(x(1-x)) dx = π
        let f = |x: f64| {
            if x > 0.0 && x < 1.0 {
                1.0 / (x * (1.0 - x)).sqrt()
            } else {
                0.0
            }
        };
        let result = integrate_tanh_sinh(&f, 0.0, 1.0, 1e-4, 200);

        assert!((result.value - PI).abs() < 0.1);
    }

    #[test]
    fn test_wynn_epsilon() {
        // Test with partial sums of alternating series: 1 - 1/2 + 1/3 - 1/4 + ... = ln(2)
        let mut partial_sums = Vec::new();
        let mut sum = 0.0;
        for n in 1..=10 {
            sum += (if n % 2 == 1 { 1.0 } else { -1.0 }) / (n as f64);
            partial_sums.push(sum);
        }

        let result = wynn_epsilon(&partial_sums);
        assert!((result - 2.0_f64.ln()).abs() < 0.01);
    }

    #[test]
    fn test_oscillatory_integral() {
        // ∫₀^{10π} sin(x)/x dx (partial Dirichlet integral)
        // Should approach π/2 as upper limit → ∞
        let result = integrate_oscillatory_real(
            &|x: f64| if x.abs() > 1e-10 { 1.0 / x } else { 1.0 },
            &|x: f64| x.sin(),
            0.001, // Avoid x=0
            10.0 * PI,
            2.0 * PI,
            1e-6,
            10,
        );

        // For finite upper bound, won't equal π/2 exactly
        assert!(result.value > 1.0);
        assert!(result.value < 2.0);
    }
}
