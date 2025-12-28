//! Polynomial integration using the power rule.
//!
//! For a polynomial p(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ,
//! the integral is:
//!
//! ∫p(x)dx = a₀x + (a₁/2)x² + (a₂/3)x³ + ... + (aₙ/(n+1))xⁿ⁺¹ + C

use tertius_poly::dense::DensePoly;
use tertius_rings::traits::Field;

/// Computes the indefinite integral of a polynomial using the power rule.
///
/// # Example
///
/// ```ignore
/// // ∫(1 + 2x + 3x²)dx = x + x² + x³
/// let p = DensePoly::new(vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)]);
/// let integral = integrate_polynomial(&p);
/// // integral = 0 + 1x + 1x² + 1x³
/// ```
pub fn integrate_polynomial<F: Field>(p: &DensePoly<F>) -> DensePoly<F> {
    if p.is_zero() {
        return DensePoly::zero();
    }

    // Result has degree one higher than input
    let mut coeffs = Vec::with_capacity(p.coeffs().len() + 1);

    // Constant term is 0 (the integration constant C)
    coeffs.push(F::zero());

    // Each coefficient aₖ becomes aₖ/(k+1) at position k+1
    for (k, coeff) in p.coeffs().iter().enumerate() {
        // Divide by (k+1)
        let divisor = F::one().mul_by_scalar((k + 1) as i64);
        let new_coeff = coeff.field_div(&divisor);
        coeffs.push(new_coeff);
    }

    DensePoly::new(coeffs)
}

/// Verifies that the integral is correct by differentiating.
///
/// Returns true if d/dx(integral) == original.
pub fn verify_polynomial_integral<F: Field>(original: &DensePoly<F>, integral: &DensePoly<F>) -> bool {
    let derivative = integral.derivative();

    // Compare coefficients
    if derivative.degree() != original.degree() {
        return false;
    }

    for i in 0..=original.degree() {
        if derivative.coeff(i) != original.coeff(i) {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn q_frac(num: i64, den: i64) -> Q {
        Q::new(num.into(), den.into())
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_integrate_constant() {
        // ∫ 5 dx = 5x
        let p = poly(&[5]);
        let integral = integrate_polynomial(&p);

        assert_eq!(integral.degree(), 1);
        assert_eq!(integral.coeff(0), q(0));
        assert_eq!(integral.coeff(1), q(5));
    }

    #[test]
    fn test_integrate_linear() {
        // ∫ (2 + 3x) dx = 2x + (3/2)x²
        let p = poly(&[2, 3]);
        let integral = integrate_polynomial(&p);

        assert_eq!(integral.degree(), 2);
        assert_eq!(integral.coeff(0), q(0));
        assert_eq!(integral.coeff(1), q(2));
        assert_eq!(integral.coeff(2), q_frac(3, 2));
    }

    #[test]
    fn test_integrate_quadratic() {
        // ∫ (1 + 2x + 3x²) dx = x + x² + x³
        let p = poly(&[1, 2, 3]);
        let integral = integrate_polynomial(&p);

        assert_eq!(integral.degree(), 3);
        assert_eq!(integral.coeff(0), q(0));
        assert_eq!(integral.coeff(1), q(1));
        assert_eq!(integral.coeff(2), q(1));
        assert_eq!(integral.coeff(3), q(1));
    }

    #[test]
    fn test_integrate_zero() {
        let p = DensePoly::<Q>::zero();
        let integral = integrate_polynomial(&p);

        assert!(integral.is_zero());
    }

    #[test]
    fn test_verify_integral() {
        let p = poly(&[1, 2, 3]);
        let integral = integrate_polynomial(&p);

        assert!(verify_polynomial_integral(&p, &integral));
    }

    #[test]
    fn test_power_rule_general() {
        // ∫ x^n dx = x^(n+1) / (n+1)
        for n in 0..5 {
            let mut coeffs = vec![q(0); n + 1];
            coeffs[n] = q(1); // x^n
            let p = DensePoly::new(coeffs);

            let integral = integrate_polynomial(&p);

            // Should be x^(n+1) / (n+1)
            assert_eq!(integral.degree(), n + 1);
            assert_eq!(integral.coeff(n + 1), q_frac(1, (n + 1) as i64));

            // Verify by differentiation
            assert!(verify_polynomial_integral(&p, &integral));
        }
    }
}
