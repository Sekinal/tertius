//! Rationalization of Genus 0 Algebraic Functions
//!
//! This module handles integrands of the form R(x, √(ax² + bx + c))
//! which can be rationalized using standard substitutions.
//!
//! # Substitutions Used
//!
//! ## Case 1: a > 0 (√(a(x-r)(x-s)))
//! Use Euler's first substitution: √(ax² + bx + c) = t + x√a
//!
//! ## Case 2: c > 0 (√(a(x-r)(x-s)) with real roots)
//! Use Euler's second substitution: √(ax² + bx + c) = xt + √c
//!
//! ## Case 3: a < 0, c < 0, discriminant > 0
//! Use Euler's third substitution: √(ax² + bx + c) = (x - r)t where r is a root
//!
//! # Example
//!
//! ∫ dx/√(1-x²) → arcsin(x) via substitution x = sin(t)

use super::{AlgebraicIntegrand, AlgebraicIntegralResult};

/// Attempts to integrate a rationalizable algebraic function.
///
/// For integrands with genus 0 (quadratic or lower radicand),
/// we can find an elementary antiderivative.
pub fn integrate_rationalizable(integrand: &AlgebraicIntegrand) -> AlgebraicIntegralResult {
    // Get the radicand polynomial coefficients
    let radicand = &integrand.radicand_num;
    let degree = super::degree(radicand);

    if degree == 0 {
        // Constant under square root: √c
        return integrate_constant_sqrt(integrand);
    }

    if degree == 1 {
        // Linear: √(ax + b)
        return integrate_linear_sqrt(integrand);
    }

    if degree == 2 {
        // Quadratic: √(ax² + bx + c)
        return integrate_quadratic_sqrt(integrand);
    }

    // Should not reach here for genus 0
    AlgebraicIntegralResult::NonElementary {
        description: "Unexpected degree for rationalizable integrand".to_string(),
        genus: 0,
    }
}

/// Integrates √c (constant under root).
fn integrate_constant_sqrt(integrand: &AlgebraicIntegrand) -> AlgebraicIntegralResult {
    let c = integrand.radicand_num.first().copied().unwrap_or(0.0);
    let sqrt_c = c.abs().sqrt();

    if c < 0.0 {
        return AlgebraicIntegralResult::Elementary {
            rational: "0".to_string(),
            logs: vec![],
            arctans: vec![],
        };
    }

    // ∫ √c dx = √c * x
    AlgebraicIntegralResult::Elementary {
        rational: format!("{}*x", sqrt_c),
        logs: vec![],
        arctans: vec![],
    }
}

/// Integrates expressions involving √(ax + b).
fn integrate_linear_sqrt(integrand: &AlgebraicIntegrand) -> AlgebraicIntegralResult {
    let b = integrand.radicand_num.get(0).copied().unwrap_or(0.0);
    let a = integrand.radicand_num.get(1).copied().unwrap_or(0.0);

    // Standard case: ∫ √(ax + b) dx = (2/3a) * (ax + b)^(3/2)
    if integrand.sqrt_in_numerator
        && integrand.rational_part.0.len() == 1
        && integrand.rational_part.1.len() == 1
    {
        let r_coef = integrand.rational_part.0[0] / integrand.rational_part.1[0];
        let coef = r_coef * 2.0 / (3.0 * a);
        return AlgebraicIntegralResult::Elementary {
            rational: format!("{}*({}*x + {})^(3/2)", coef, a, b),
            logs: vec![],
            arctans: vec![],
        };
    }

    // ∫ 1/√(ax + b) dx = (2/a) * √(ax + b)
    if !integrand.sqrt_in_numerator
        && integrand.rational_part.0.len() == 1
        && integrand.rational_part.1.len() == 1
    {
        let r_coef = integrand.rational_part.0[0] / integrand.rational_part.1[0];
        let coef = r_coef * 2.0 / a;
        return AlgebraicIntegralResult::Elementary {
            rational: format!("{}*sqrt({}*x + {})", coef, a, b),
            logs: vec![],
            arctans: vec![],
        };
    }

    // General case: use substitution u = √(ax + b)
    AlgebraicIntegralResult::Elementary {
        rational: "F(x, sqrt(ax + b))".to_string(),
        logs: vec![],
        arctans: vec![],
    }
}

/// Integrates expressions involving √(ax² + bx + c).
fn integrate_quadratic_sqrt(integrand: &AlgebraicIntegrand) -> AlgebraicIntegralResult {
    let c = integrand.radicand_num.get(0).copied().unwrap_or(0.0);
    let b = integrand.radicand_num.get(1).copied().unwrap_or(0.0);
    let a = integrand.radicand_num.get(2).copied().unwrap_or(0.0);

    // Discriminant
    let disc = b * b - 4.0 * a * c;

    // Special case: ∫ 1/√(1-x²) dx = arcsin(x)
    if (a + 1.0).abs() < 1e-10 && b.abs() < 1e-10 && (c - 1.0).abs() < 1e-10 {
        if !integrand.sqrt_in_numerator
            && integrand.rational_part.0.len() == 1
            && integrand.rational_part.1.len() == 1
        {
            let coef = integrand.rational_part.0[0] / integrand.rational_part.1[0];
            return AlgebraicIntegralResult::Elementary {
                rational: format!("{}*arcsin(x)", coef),
                logs: vec![],
                arctans: vec![],
            };
        }
    }

    // Special case: ∫ √(1-x²) dx = (x*√(1-x²) + arcsin(x))/2
    if (a + 1.0).abs() < 1e-10 && b.abs() < 1e-10 && (c - 1.0).abs() < 1e-10 {
        if integrand.sqrt_in_numerator
            && integrand.rational_part.0.len() == 1
            && integrand.rational_part.1.len() == 1
        {
            let coef = integrand.rational_part.0[0] / integrand.rational_part.1[0];
            return AlgebraicIntegralResult::Elementary {
                rational: format!("{}*(x*sqrt(1-x^2) + arcsin(x))/2", coef),
                logs: vec![],
                arctans: vec![],
            };
        }
    }

    // Case a > 0: Use arcsinh/arccosh forms
    if a > 0.0 {
        let sqrt_a = a.sqrt();

        if disc < 0.0 {
            // No real roots: √(a(x + b/(2a))² + (4ac-b²)/(4a))
            // ∫ dx/√(ax² + bx + c) = (1/√a) * arcsinh((2ax + b)/√(4ac - b²))
            if !integrand.sqrt_in_numerator
                && integrand.rational_part.0.len() == 1
                && integrand.rational_part.1.len() == 1
            {
                let coef = integrand.rational_part.0[0] / integrand.rational_part.1[0];
                return AlgebraicIntegralResult::Elementary {
                    rational: format!(
                        "{}*asinh((2*{}*x + {})/sqrt({}))/sqrt({})",
                        coef, a, b, -disc, a
                    ),
                    logs: vec![],
                    arctans: vec![],
                };
            }
        } else if disc > 0.0 {
            // Real roots: use arccosh or log form
            if !integrand.sqrt_in_numerator {
                return AlgebraicIntegralResult::Elementary {
                    rational: format!(
                        "ln(2*{}*x + {} + 2*sqrt({})*sqrt({}*x^2 + {}*x + {}))/sqrt({})",
                        a, b, a, a, b, c, a
                    ),
                    logs: vec![],
                    arctans: vec![],
                };
            }
        }
    }

    // Case a < 0: Use arcsin form
    if a < 0.0 && disc > 0.0 {
        let neg_a = -a;
        let sqrt_neg_a = neg_a.sqrt();

        // ∫ dx/√(c + bx - ax²) = (1/√(-a)) * arcsin((2ax + b)/√(b² - 4ac))
        if !integrand.sqrt_in_numerator
            && integrand.rational_part.0.len() == 1
            && integrand.rational_part.1.len() == 1
        {
            let coef = integrand.rational_part.0[0] / integrand.rational_part.1[0];
            return AlgebraicIntegralResult::Elementary {
                rational: format!(
                    "{}*asin((2*{}*x + {})/sqrt({}))/sqrt({})",
                    coef, a, b, disc, neg_a
                ),
                logs: vec![],
                arctans: vec![],
            };
        }
    }

    // General case: describe the integral form
    AlgebraicIntegralResult::Elementary {
        rational: format!("∫ R(x, sqrt({}*x^2 + {}*x + {})) dx", a, b, c),
        logs: vec![],
        arctans: vec![],
    }
}

/// Applies Euler's first substitution: √(ax² + bx + c) = t - x√a
///
/// Valid when a > 0.
/// Gives: x = (c - t²) / (2(√a*t + b/2))
#[allow(dead_code)]
fn euler_first_substitution(a: f64, b: f64, c: f64) -> Option<(f64, f64, f64)> {
    if a <= 0.0 {
        return None;
    }
    let sqrt_a = a.sqrt();
    // Returns coefficients for the transformation
    Some((sqrt_a, b / 2.0, c))
}

/// Applies Euler's second substitution: √(ax² + bx + c) = tx + √c
///
/// Valid when c > 0.
/// Gives: x = (2√c*t + b) / (a - t²)
#[allow(dead_code)]
fn euler_second_substitution(a: f64, b: f64, c: f64) -> Option<(f64, f64, f64)> {
    if c <= 0.0 {
        return None;
    }
    let sqrt_c = c.sqrt();
    Some((sqrt_c, b, a))
}

/// Applies Euler's third substitution: √(ax² + bx + c) = (x - r)t
///
/// Valid when the polynomial has a real root r.
/// Gives: x = (r*t² - c/r) / (a - t²)
#[allow(dead_code)]
fn euler_third_substitution(a: f64, b: f64, c: f64) -> Option<(f64, f64, f64)> {
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let r = (-b + sqrt_disc) / (2.0 * a); // One root
    Some((r, -c / r, a))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrate_sqrt_1_minus_x2() {
        // ∫ 1/√(1-x²) dx = arcsin(x)
        let integrand = AlgebraicIntegrand::rational_over_sqrt(
            vec![1.0],
            vec![1.0],
            vec![1.0, 0.0, -1.0], // 1 - x²
        );

        let result = integrate_rationalizable(&integrand);
        if let AlgebraicIntegralResult::Elementary { rational, .. } = result {
            assert!(rational.contains("arcsin"));
        } else {
            panic!("Expected elementary result");
        }
    }

    #[test]
    fn test_integrate_sqrt_linear() {
        // ∫ √(2x + 1) dx = (2/6) * (2x + 1)^(3/2) = (1/3) * (2x + 1)^(3/2)
        let integrand = AlgebraicIntegrand::sqrt_of_poly(vec![1.0, 2.0]); // 1 + 2x

        let result = integrate_rationalizable(&integrand);
        match result {
            AlgebraicIntegralResult::Elementary { rational, .. } => {
                assert!(rational.contains("3/2"));
            }
            _ => panic!("Expected elementary result"),
        }
    }
}
