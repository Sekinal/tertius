//! Integration of √(P(x)/Q(x)) Type Functions
//!
//! This module handles integrands involving square roots of rational functions,
//! particularly the form √((1+x)/(1-x)) which appears in the golden ratio integral.
//!
//! # Key Substitutions
//!
//! For √((1+x)/(1-x)):
//! - Let t = √((1+x)/(1-x)), then t² = (1+x)/(1-x)
//! - Solving: x = (t² - 1)/(t² + 1)
//! - dx = 4t/(t² + 1)² dt
//!
//! This transforms the integral into a rational function of t.

use super::AlgebraicIntegralResult;

/// Transforms √((1+x)/(1-x)) using the substitution t = √((1+x)/(1-x)).
///
/// Returns the transformed integrand as a rational function of t.
///
/// x = (t² - 1)/(t² + 1)
/// dx = 4t/(t² + 1)² dt
/// √((1+x)/(1-x)) = t
pub fn sqrt_ratio_substitution(x: f64) -> (f64, f64) {
    // t = √((1+x)/(1-x))
    let t = ((1.0 + x) / (1.0 - x)).sqrt();
    // dx/dt = 4t/(t² + 1)²
    let dt_dx = (t * t + 1.0).powi(2) / (4.0 * t);
    (t, dt_dx)
}

/// Inverse of the sqrt ratio substitution.
///
/// Given t, returns x = (t² - 1)/(t² + 1).
pub fn sqrt_ratio_inverse(t: f64) -> f64 {
    let t2 = t * t;
    (t2 - 1.0) / (t2 + 1.0)
}

/// For the golden ratio integral, after substitution t = √((1+x)/(1-x)),
/// the integral becomes:
///
/// ∫ (1/x) * t * ln((2x²+2x+1)/(2x²-2x+1)) dx
///
/// With x = (t²-1)/(t²+1), we substitute and get a rational integrand in t.
pub fn transform_golden_ratio_integrand(t: f64) -> f64 {
    // x in terms of t
    let t2 = t * t;
    let x = (t2 - 1.0) / (t2 + 1.0);

    if x.abs() < 1e-10 {
        // Limit at x=0 (t=1)
        return 4.0 * 4.0 / (t2 + 1.0).powi(2);
    }

    // The logarithm term
    let ln_num = 2.0 * x * x + 2.0 * x + 1.0;
    let ln_den = 2.0 * x * x - 2.0 * x + 1.0;
    let ln_term = (ln_num / ln_den).ln();

    // (1/x) * t * ln(...) * (4t/(t²+1)²)
    // = (4t²/(x*(t²+1)²)) * ln(...)
    (4.0 * t2 / (x * (t2 + 1.0).powi(2))) * ln_term
}

/// Represents the parameters for a Möbius transformation.
///
/// Transforms x → (ax + b)/(cx + d)
#[derive(Clone, Debug)]
pub struct MobiusTransform {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl MobiusTransform {
    /// Creates a new Möbius transformation.
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }

    /// Applies the transformation to x.
    pub fn apply(&self, x: f64) -> f64 {
        (self.a * x + self.b) / (self.c * x + self.d)
    }

    /// Computes the derivative of the transformation.
    ///
    /// d/dx[(ax+b)/(cx+d)] = (ad - bc)/(cx + d)²
    pub fn derivative(&self, x: f64) -> f64 {
        let det = self.a * self.d - self.b * self.c;
        let denom = self.c * x + self.d;
        det / (denom * denom)
    }

    /// Computes the inverse transformation.
    ///
    /// If T(x) = (ax+b)/(cx+d), then T⁻¹(y) = (dy - b)/(-cy + a)
    pub fn inverse(&self) -> Self {
        Self {
            a: self.d,
            b: -self.b,
            c: -self.c,
            d: self.a,
        }
    }
}

/// For integrands of the form (1/x) * √((1+x)/(1-x)) * f(x),
/// apply the substitution and return integration result.
pub fn integrate_sqrt_ratio_times_func<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
) -> AlgebraicIntegralResult {
    // Transform the integral using t = √((1+x)/(1-x))
    // The bounds transform as:
    // x = a → t_a = √((1+a)/(1-a))
    // x = b → t_b = √((1+b)/(1-b))

    // Handle endpoint singularities
    let a_adjusted = if a <= -1.0 { -1.0 + 1e-10 } else { a };
    let b_adjusted = if b >= 1.0 { 1.0 - 1e-10 } else { b };

    // Transform bounds
    let t_a = if a_adjusted >= 1.0 {
        f64::INFINITY
    } else {
        ((1.0 + a_adjusted) / (1.0 - a_adjusted)).sqrt()
    };

    let t_b = if b_adjusted >= 1.0 {
        f64::INFINITY
    } else {
        ((1.0 + b_adjusted) / (1.0 - b_adjusted)).sqrt()
    };

    // The transformed integrand in t
    let transformed = |t: f64| {
        let t2 = t * t;
        let x = (t2 - 1.0) / (t2 + 1.0);
        let jacobian = 4.0 * t / (t2 + 1.0).powi(2);

        if x.abs() < 1e-12 {
            0.0 // Handle x=0 specially
        } else {
            (t / x) * f(x) * jacobian
        }
    };

    // Integrate in t-space
    use crate::numerical::adaptive_integrate;

    let result = if t_a.is_finite() && t_b.is_finite() {
        adaptive_integrate(&transformed, t_a, t_b, tolerance, tolerance, 1000)
    } else {
        // Handle infinite bounds using another substitution
        use crate::numerical::adaptive_integrate_to_infinity;
        if !t_a.is_finite() && !t_b.is_finite() {
            use crate::numerical::adaptive_integrate_full_line;
            adaptive_integrate_full_line(&transformed, tolerance, tolerance, 1000)
        } else if !t_b.is_finite() {
            adaptive_integrate_to_infinity(&transformed, t_a.max(0.0), tolerance, tolerance, 1000)
        } else {
            // t_a is infinite
            use crate::numerical::adaptive_integrate_from_neg_infinity;
            adaptive_integrate_from_neg_infinity(&transformed, t_b, tolerance, tolerance, 1000)
        }
    };

    AlgebraicIntegralResult::Numerical(result.value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sqrt_ratio_substitution() {
        // At x = 0: t = 1
        let (t, _) = sqrt_ratio_substitution(0.0);
        assert!((t - 1.0).abs() < 1e-10);

        // At x = 0.5: t = √3
        let (t, _) = sqrt_ratio_substitution(0.5);
        assert!((t - 3.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_sqrt_ratio_inverse() {
        // t = 1 → x = 0
        let x = sqrt_ratio_inverse(1.0);
        assert!(x.abs() < 1e-10);

        // t = √3 → x = 0.5
        let x = sqrt_ratio_inverse(3.0_f64.sqrt());
        assert!((x - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_mobius_transform() {
        // Identity transform: x → x
        let id = MobiusTransform::new(1.0, 0.0, 0.0, 1.0);
        assert!((id.apply(2.0) - 2.0).abs() < 1e-10);

        // x → (x-1)/(x+1), inverse of t² = (1+x)/(1-x) up to squares
        let t = MobiusTransform::new(1.0, -1.0, 1.0, 1.0);
        // At x = 0: (0-1)/(0+1) = -1
        assert!((t.apply(0.0) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_transform_at_origin() {
        // At x = 0 (t = 1), the transformed integrand should be finite
        let val = transform_golden_ratio_integrand(1.0);
        assert!(val.is_finite());
    }
}
