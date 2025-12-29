//! Beta and Gamma Function Integrals
//!
//! This module handles integrals related to the Beta and Gamma functions:
//!
//! - ∫₀^1 x^(a-1)(1-x)^(b-1) dx = B(a,b) = Γ(a)Γ(b)/Γ(a+b)
//! - ∫₀^∞ x^(n-1)e^(-x) dx = Γ(n)
//! - ∫₀^∞ x^(a-1)e^(-bx) dx = Γ(a)/b^a

use std::f64::consts::PI;
use super::{gamma, beta};

/// Beta function integral evaluator.
#[derive(Clone, Debug)]
pub struct BetaIntegral {
    /// First parameter a
    pub a: f64,
    /// Second parameter b
    pub b: f64,
}

impl BetaIntegral {
    /// Creates a new Beta integral B(a, b).
    pub fn new(a: f64, b: f64) -> Self {
        Self { a, b }
    }

    /// Evaluates the Beta function B(a, b) = ∫₀^1 x^(a-1)(1-x)^(b-1) dx.
    pub fn evaluate(&self) -> f64 {
        beta(self.a, self.b)
    }

    /// The integrand x^(a-1)(1-x)^(b-1).
    pub fn integrand(&self, x: f64) -> f64 {
        if x <= 0.0 || x >= 1.0 {
            return 0.0;
        }
        x.powf(self.a - 1.0) * (1.0 - x).powf(self.b - 1.0)
    }

    /// Evaluates numerically for verification.
    pub fn evaluate_numerical(&self, tolerance: f64) -> f64 {
        use crate::numerical::adaptive_integrate;

        // Handle endpoint singularities
        let epsilon = 1e-10;
        let result = adaptive_integrate(
            &|x| self.integrand(x),
            epsilon,
            1.0 - epsilon,
            tolerance,
            tolerance,
            1000,
        );
        result.value
    }
}

/// Gamma function integral evaluator.
#[derive(Clone, Debug)]
pub struct GammaIntegral {
    /// Parameter n
    pub n: f64,
}

impl GammaIntegral {
    /// Creates a new Gamma integral Γ(n).
    pub fn new(n: f64) -> Self {
        Self { n }
    }

    /// Evaluates Γ(n) = ∫₀^∞ x^(n-1)e^(-x) dx.
    pub fn evaluate(&self) -> f64 {
        gamma(self.n)
    }

    /// The integrand x^(n-1)e^(-x).
    pub fn integrand(&self, x: f64) -> f64 {
        if x <= 0.0 {
            if self.n >= 1.0 {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        }
        x.powf(self.n - 1.0) * (-x).exp()
    }

    /// Evaluates numerically for verification.
    pub fn evaluate_numerical(&self, tolerance: f64) -> f64 {
        use crate::numerical::adaptive_integrate_to_infinity;

        let epsilon = if self.n >= 1.0 { 0.0 } else { 1e-10 };
        let result = adaptive_integrate_to_infinity(
            &|x| self.integrand(x),
            epsilon,
            tolerance,
            tolerance,
            1000,
        );
        result.value
    }
}

/// Recognizes and evaluates common Beta/Gamma integral patterns.
pub fn recognize_beta_pattern<F: Fn(f64) -> f64>(
    f: &F,
    a_bound: f64,
    b_bound: f64,
) -> Option<BetaIntegral> {
    // Check if bounds are [0, 1]
    if (a_bound - 0.0).abs() > 1e-10 || (b_bound - 1.0).abs() > 1e-10 {
        return None;
    }

    // Sample the function to estimate parameters
    // For x^(a-1)(1-x)^(b-1), taking logs:
    // ln(f(x)) = (a-1)ln(x) + (b-1)ln(1-x)

    // Sample at x = 0.25 and x = 0.75
    let f1 = f(0.25);
    let f2 = f(0.75);

    if f1 <= 0.0 || f2 <= 0.0 {
        return None;
    }

    let ln_f1 = f1.ln();
    let ln_f2 = f2.ln();

    // At x = 0.25: ln(f) = (a-1)ln(0.25) + (b-1)ln(0.75)
    // At x = 0.75: ln(f) = (a-1)ln(0.75) + (b-1)ln(0.25)

    let ln_quarter: f64 = 0.25_f64.ln();
    let ln_three_quarter: f64 = 0.75_f64.ln();

    // Solve the system:
    // (a-1)*ln(0.25) + (b-1)*ln(0.75) = ln_f1
    // (a-1)*ln(0.75) + (b-1)*ln(0.25) = ln_f2

    let det = ln_quarter * ln_quarter - ln_three_quarter * ln_three_quarter;
    if det.abs() < 1e-10 {
        return None;
    }

    let a_minus_1 = (ln_f1 * ln_quarter - ln_f2 * ln_three_quarter) / det;
    let b_minus_1 = (ln_f2 * ln_quarter - ln_f1 * ln_three_quarter) / det;

    let a = a_minus_1 + 1.0;
    let b = b_minus_1 + 1.0;

    // Verify at another point
    let test_x: f64 = 0.5;
    let expected = test_x.powf(a - 1.0) * (1.0 - test_x).powf(b - 1.0);
    let actual = f(test_x);

    if (expected - actual).abs() / expected.abs().max(1.0) < 0.01 {
        Some(BetaIntegral::new(a, b))
    } else {
        None
    }
}

/// Wallis integral: ∫₀^(π/2) sin^n(x) dx.
///
/// Related to Beta function: = (1/2)B((n+1)/2, 1/2)
pub fn wallis_integral(n: u32) -> f64 {
    let a = (n as f64 + 1.0) / 2.0;
    let b = 0.5;
    0.5 * beta(a, b)
}

/// Computes factorials using Gamma function.
pub fn factorial(n: u64) -> f64 {
    gamma(n as f64 + 1.0)
}

/// Computes the double factorial n!! using Gamma function.
pub fn double_factorial(n: u64) -> f64 {
    if n == 0 || n == 1 {
        return 1.0;
    }

    let k = n as f64;
    if n % 2 == 0 {
        // n!! = 2^(n/2) * (n/2)!
        2.0_f64.powf(k / 2.0) * gamma(k / 2.0 + 1.0)
    } else {
        // n!! = n! / (n-1)!!
        // Using: n!! = 2^((n+1)/2) * Γ((n+2)/2) / √π
        2.0_f64.powf((k + 1.0) / 2.0) * gamma((k + 2.0) / 2.0) / PI.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_beta_integral() {
        let b = BetaIntegral::new(2.0, 3.0);
        let exact = b.evaluate();
        let numerical = b.evaluate_numerical(1e-6);

        assert!((exact - numerical).abs() < 1e-4);
        assert!((exact - 1.0 / 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_gamma_integral() {
        let g = GammaIntegral::new(4.0);
        let exact = g.evaluate();

        assert!((exact - 6.0).abs() < 1e-10); // 3!
    }

    #[test]
    fn test_wallis_integral() {
        // ∫₀^(π/2) sin²(x) dx = π/4
        let result = wallis_integral(2);
        assert!((result - PI / 4.0).abs() < 1e-10);

        // ∫₀^(π/2) sin³(x) dx = 2/3
        let result3 = wallis_integral(3);
        assert!((result3 - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_factorial() {
        assert!((factorial(0) - 1.0).abs() < 1e-10);
        assert!((factorial(1) - 1.0).abs() < 1e-10);
        assert!((factorial(5) - 120.0).abs() < 1e-10);
    }

    #[test]
    fn test_double_factorial() {
        // 5!! = 5 * 3 * 1 = 15
        assert!((double_factorial(5) - 15.0).abs() < 1e-10);

        // 6!! = 6 * 4 * 2 = 48
        assert!((double_factorial(6) - 48.0).abs() < 1e-10);
    }

    #[test]
    fn test_recognize_beta() {
        // f(x) = x(1-x) = x^1 * (1-x)^1, so a=2, b=2
        let f = |x: f64| x * (1.0 - x);
        let result = recognize_beta_pattern(&f, 0.0, 1.0);

        assert!(result.is_some());
        let beta = result.unwrap();
        assert!((beta.a - 2.0).abs() < 0.1);
        assert!((beta.b - 2.0).abs() < 0.1);
    }
}
