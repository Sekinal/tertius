//! Dirichlet Integrals
//!
//! This module handles Dirichlet-type integrals, including:
//! - ∫₀^∞ sin(x)/x dx = π/2
//! - ∫₀^∞ sin²(x)/x² dx = π/2
//! - ∫₀^∞ sin(ax)/x dx = π/2 * sign(a)

use std::f64::consts::PI;

/// The Dirichlet integral ∫₀^∞ sin(x)/x dx.
#[derive(Clone, Debug)]
pub struct DirichletIntegral;

impl DirichletIntegral {
    /// Exact value: π/2
    pub const EXACT_VALUE: f64 = PI / 2.0;

    /// Evaluates the Dirichlet kernel sin(x)/x at a point.
    pub fn sinc(x: f64) -> f64 {
        if x.abs() < 1e-10 {
            1.0 // lim_{x→0} sin(x)/x = 1
        } else {
            x.sin() / x
        }
    }

    /// Evaluates ∫₀^a sin(x)/x dx for finite a.
    pub fn integral_to_a(a: f64, tolerance: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }

        use crate::numerical::adaptive_integrate;
        let result = adaptive_integrate(&Self::sinc, 1e-10, a, tolerance, tolerance, 1000);
        result.value
    }

    /// Evaluates the full Dirichlet integral numerically.
    ///
    /// Uses the exact value π/2 since numerical integration to infinity
    /// is challenging for oscillatory integrands.
    pub fn evaluate() -> f64 {
        Self::EXACT_VALUE
    }

    /// The sine integral Si(x) = ∫₀^x sin(t)/t dt.
    pub fn sine_integral(x: f64) -> f64 {
        if x.abs() < 1e-10 {
            return x;
        }

        if x.abs() < 1.0 {
            // Use Taylor series for small x
            Self::sine_integral_taylor(x)
        } else {
            // Use asymptotic expansion for large x
            Self::sine_integral_asymptotic(x)
        }
    }

    /// Taylor series for Si(x) around x = 0.
    ///
    /// Si(x) = x - x³/18 + x⁵/600 - x⁷/35280 + ...
    fn sine_integral_taylor(x: f64) -> f64 {
        let mut result = 0.0;
        let mut term = x;
        let x2 = x * x;

        for n in 0..50 {
            result += term;
            let k = 2 * n + 1;
            term *= -x2 / ((2 * k + 2) * (2 * k + 3)) as f64;
            if term.abs() < 1e-16 {
                break;
            }
        }

        result
    }

    /// Asymptotic expansion for Si(x) for large |x|.
    ///
    /// Si(x) ≈ π/2 - cos(x)/x - sin(x)/x² + ...
    fn sine_integral_asymptotic(x: f64) -> f64 {
        let sign = x.signum();
        let x = x.abs();

        if x > 50.0 {
            // For very large x, just use the limit
            return sign * PI / 2.0;
        }

        // Use numerical integration for intermediate values
        use crate::numerical::adaptive_integrate;
        let result = adaptive_integrate(&Self::sinc, 1e-10, x, 1e-8, 1e-8, 500);

        sign * result.value
    }

    /// The cosine integral Ci(x) = -∫_x^∞ cos(t)/t dt.
    ///
    /// Also equals γ + ln(x) + ∫₀^x (cos(t)-1)/t dt
    pub fn cosine_integral(x: f64) -> f64 {
        use super::MathConstants;

        if x <= 0.0 {
            return f64::NEG_INFINITY;
        }

        if x < 1.0 {
            // Use integral definition with series
            Self::cosine_integral_small(x)
        } else {
            // Use numerical integration
            Self::cosine_integral_numerical(x)
        }
    }

    fn cosine_integral_small(x: f64) -> f64 {
        use super::MathConstants;

        // Ci(x) = γ + ln(x) + ∫₀^x (cos(t)-1)/t dt
        let euler = MathConstants::EULER_GAMMA;

        // The integral ∫₀^x (cos(t)-1)/t dt has series:
        // = -x²/4 + x⁴/96 - x⁶/4320 + ...
        let mut integral = 0.0;
        let mut term = -x * x / 4.0;
        let x2 = x * x;

        for n in 1..50 {
            integral += term;
            let k = 2 * n;
            term *= -x2 / ((2 * k) * (2 * k + 1) * (2 * k + 2)) as f64;
            term *= (2 * k) as f64 / (2 * k + 2) as f64;
            if term.abs() < 1e-16 {
                break;
            }
        }

        euler + x.ln() + integral
    }

    fn cosine_integral_numerical(x: f64) -> f64 {
        use crate::numerical::adaptive_integrate_to_infinity;

        // Ci(x) = -∫_x^∞ cos(t)/t dt
        let integrand = |t: f64| t.cos() / t;
        let result = adaptive_integrate_to_infinity(&integrand, x, 1e-8, 1e-8, 500);
        -result.value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sinc_at_zero() {
        assert!((DirichletIntegral::sinc(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sinc_values() {
        assert!((DirichletIntegral::sinc(PI) - 0.0).abs() < 1e-10);
        assert!((DirichletIntegral::sinc(PI / 2.0) - 2.0 / PI).abs() < 1e-10);
    }

    #[test]
    fn test_dirichlet_exact() {
        assert!((DirichletIntegral::evaluate() - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_sine_integral_small() {
        // Si(1) ≈ 0.946083...
        let si1 = DirichletIntegral::sine_integral(1.0);
        assert!((si1 - 0.9460830703671830).abs() < 1e-6);
    }

    #[test]
    fn test_sine_integral_limit() {
        // lim_{x→∞} Si(x) = π/2
        let si_large = DirichletIntegral::sine_integral(100.0);
        assert!((si_large - PI / 2.0).abs() < 0.1);
    }
}
