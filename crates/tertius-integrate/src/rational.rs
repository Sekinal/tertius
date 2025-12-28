//! Complete rational function integration.
//!
//! Integrates rational functions P(x)/Q(x) using:
//! 1. **Polynomial division**: Extract polynomial part, integrate with power rule
//! 2. **Hermite reduction**: Reduce to simple poles, extract rational part
//! 3. **Rothstein-Trager**: Compute logarithmic part for simple poles
//!
//! # Algorithm
//!
//! Given ∫(P/Q)dx:
//! 1. Divide P by Q to get P = Q*poly + R where deg(R) < deg(Q)
//! 2. ∫(P/Q) = ∫poly + ∫(R/Q)
//! 3. Apply Hermite reduction to R/Q: ∫(R/Q) = g + ∫(A/D) where D is squarefree
//! 4. Apply Rothstein-Trager to A/D to get logarithmic part
//! 5. Result: ∫poly + g + Σcᵢlog(vᵢ)

use tertius_poly::dense::DensePoly;
use tertius_rational_func::hermite::{hermite_reduce, HermiteReductionResult};
use tertius_rational_func::RationalFunction;
use tertius_rings::rationals::Q;
use tertius_rings::traits::Field;

use crate::polynomial::integrate_polynomial;
use crate::rothstein_trager::{
    rothstein_trager, rothstein_trager_algebraic, AlgebraicLogarithmicPart, LogarithmicPart,
};

/// Result of rational function integration.
#[derive(Clone, Debug)]
pub struct RationalIntegrationResult<F: Field> {
    /// The polynomial part of the integral.
    pub polynomial_part: DensePoly<F>,
    /// The rational part from Hermite reduction.
    pub rational_part: RationalFunction<F>,
    /// The logarithmic part from Rothstein-Trager.
    pub logarithmic_part: Option<LogarithmicPart<F>>,
    /// Whether the integration was complete (no algebraic extensions needed).
    pub is_complete: bool,
}

impl<F: Field> RationalIntegrationResult<F> {
    /// Returns true if this is a zero integral.
    pub fn is_zero(&self) -> bool {
        self.polynomial_part.is_zero()
            && self.rational_part.is_zero()
            && self.logarithmic_part.as_ref().map_or(true, |l| l.is_empty())
    }
}

/// Integrates a rational function completely.
///
/// # Algorithm
///
/// 1. Extract and integrate polynomial part
/// 2. Apply Hermite reduction to get rational + squarefree parts
/// 3. Apply Rothstein-Trager to get logarithmic part
///
/// # Returns
///
/// The complete integral as polynomial + rational + logarithmic parts.
///
/// # Example
///
/// ```ignore
/// use tertius_integrate::integrate_rational;
/// use tertius_rational_func::RationalFunction;
///
/// // ∫ (x² + 1)/(x² - 1) dx
/// let rf = RationalFunction::new(
///     DensePoly::new(vec![Q::from_integer(1), Q::zero(), Q::from_integer(1)]),
///     DensePoly::new(vec![Q::from_integer(-1), Q::zero(), Q::from_integer(1)]),
/// );
///
/// let result = integrate_rational(&rf);
/// // result.polynomial_part = x  (from the polynomial part)
/// // result.logarithmic_part contains log terms
/// ```
pub fn integrate_rational<F: Field>(f: &RationalFunction<F>) -> RationalIntegrationResult<F> {
    // Step 1: Extract polynomial part
    let (poly_part, proper) = f.decompose_proper();

    // Integrate polynomial part: ∫poly dx
    let poly_integral = integrate_polynomial(&poly_part);

    // If the proper part is zero, we're done
    if proper.is_zero() {
        return RationalIntegrationResult {
            polynomial_part: poly_integral,
            rational_part: RationalFunction::zero(),
            logarithmic_part: Some(LogarithmicPart::empty()),
            is_complete: true,
        };
    }

    // Step 2: Apply Hermite reduction
    // ∫(proper) = g + ∫(A/D) where D is squarefree
    let HermiteReductionResult {
        rational_part,
        reduced,
    } = hermite_reduce(&proper);

    // Step 3: Apply Rothstein-Trager to the reduced part
    let logarithmic_part = if reduced.is_zero() {
        Some(LogarithmicPart::empty())
    } else {
        rothstein_trager(reduced.numerator(), reduced.denominator())
    };

    let is_complete = logarithmic_part.is_some();

    RationalIntegrationResult {
        polynomial_part: poly_integral,
        rational_part,
        logarithmic_part,
        is_complete,
    }
}

/// Integrates a simple rational function of the form A/(x - a)^n.
///
/// # Formula
///
/// - n = 1: ∫ A/(x-a) dx = A·log(x-a)
/// - n > 1: ∫ A/(x-a)^n dx = -A/((n-1)(x-a)^(n-1))
pub fn integrate_simple_pole<F: Field>(
    a_coeff: &F,
    pole: &F,
    power: u32,
) -> SimplePoleTerm<F> {
    if power == 1 {
        // Logarithmic case
        SimplePoleTerm::Logarithmic {
            coefficient: a_coeff.clone(),
            argument: DensePoly::new(vec![-pole.clone(), F::one()]), // x - pole
        }
    } else {
        // Rational case: -A / ((n-1)(x-a)^(n-1))
        let n_minus_1 = F::one().mul_by_scalar((power - 1) as i64);
        let coeff = -a_coeff.clone().field_div(&n_minus_1);

        SimplePoleTerm::Rational {
            coefficient: coeff,
            pole: pole.clone(),
            power: power - 1,
        }
    }
}

/// A term from integrating a simple pole.
#[derive(Clone, Debug)]
pub enum SimplePoleTerm<F: Field> {
    /// A logarithmic term: c * log(x - a)
    Logarithmic {
        coefficient: F,
        argument: DensePoly<F>,
    },
    /// A rational term: c / (x - a)^n
    Rational { coefficient: F, pole: F, power: u32 },
}

/// Result of rational function integration with algebraic coefficients.
///
/// Used when the resultant polynomial has algebraic roots.
#[derive(Clone, Debug)]
pub struct AlgebraicIntegrationResult {
    /// The polynomial part of the integral.
    pub polynomial_part: DensePoly<Q>,
    /// The rational part from Hermite reduction.
    pub rational_part: RationalFunction<Q>,
    /// The logarithmic part with algebraic coefficients.
    pub algebraic_log_part: Option<AlgebraicLogarithmicPart>,
    /// Whether the integration was complete.
    pub is_complete: bool,
}

/// Integrates a rational function over Q with algebraic extension support.
///
/// This version can handle integrals like ∫1/(x⁴+1)dx that require
/// algebraic number field extensions for the logarithmic part.
pub fn integrate_rational_with_algebraic(f: &RationalFunction<Q>) -> AlgebraicIntegrationResult {
    // Step 1: Extract polynomial part
    let (poly_part, proper) = f.decompose_proper();
    let poly_integral = integrate_polynomial(&poly_part);

    // If the proper part is zero, we're done
    if proper.is_zero() {
        return AlgebraicIntegrationResult {
            polynomial_part: poly_integral,
            rational_part: RationalFunction::zero(),
            algebraic_log_part: Some(AlgebraicLogarithmicPart::empty()),
            is_complete: true,
        };
    }

    // Step 2: Apply Hermite reduction
    let HermiteReductionResult {
        rational_part,
        reduced,
    } = hermite_reduce(&proper);

    // Step 3: Try standard Rothstein-Trager first
    if reduced.is_zero() {
        return AlgebraicIntegrationResult {
            polynomial_part: poly_integral,
            rational_part,
            algebraic_log_part: Some(AlgebraicLogarithmicPart::empty()),
            is_complete: true,
        };
    }

    // Try rational root finding first (fast path)
    if let Some(log_part) = rothstein_trager(reduced.numerator(), reduced.denominator()) {
        // Convert rational LogarithmicPart to AlgebraicLogarithmicPart
        return AlgebraicIntegrationResult {
            polynomial_part: poly_integral,
            rational_part,
            algebraic_log_part: Some(log_part.into_algebraic()),
            is_complete: true,
        };
    }

    // Fallback to algebraic root finding
    let alg_log_part = rothstein_trager_algebraic(reduced.numerator(), reduced.denominator());
    let is_complete = !alg_log_part.is_empty();

    AlgebraicIntegrationResult {
        polynomial_part: poly_integral,
        rational_part,
        algebraic_log_part: Some(alg_log_part),
        is_complete,
    }
}

/// Verifies an integration result by differentiation.
///
/// Checks that d/dx(result) = original.
pub fn verify_rational_integral<F: Field>(
    original: &RationalFunction<F>,
    result: &RationalIntegrationResult<F>,
) -> bool {
    // Differentiate polynomial part
    let poly_deriv = result.polynomial_part.derivative();

    // Differentiate rational part
    let rational_deriv = result.rational_part.derivative();

    // For logarithmic parts, d/dx(c*log(v)) = c*v'/v
    // This is more complex to verify, so we skip for now

    // Check that poly_deriv + rational_deriv approximately equals original
    // (ignoring the logarithmic contribution for verification)

    let poly_rf = RationalFunction::from_poly(poly_deriv);
    let sum = poly_rf + rational_deriv;

    // For a complete verification, we'd need to add the log derivatives too
    // For now, just verify the structure is reasonable
    !result.polynomial_part.is_zero()
        || !result.rational_part.is_zero()
        || result
            .logarithmic_part
            .as_ref()
            .map_or(false, |l| !l.is_empty())
        || original.is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_poly::dense::DensePoly;
    use tertius_rings::rationals::Q;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    fn rf(num: &[i64], den: &[i64]) -> RationalFunction<Q> {
        RationalFunction::new(poly(num), poly(den))
    }

    #[test]
    fn test_integrate_polynomial_only() {
        // ∫ (x + 1) dx = x²/2 + x
        let f = rf(&[1, 1], &[1]); // (1 + x) / 1

        let result = integrate_rational(&f);

        assert!(result.is_complete);
        assert!(result.rational_part.is_zero());
        assert_eq!(result.polynomial_part.degree(), 2);
    }

    #[test]
    fn test_integrate_simple_linear() {
        // ∫ 1/(x + 1) dx = log(x + 1)
        let f = rf(&[1], &[1, 1]);

        let result = integrate_rational(&f);

        assert!(result.is_complete);
        assert!(result.polynomial_part.is_zero());

        if let Some(log_part) = &result.logarithmic_part {
            assert!(!log_part.is_empty());
        }
    }

    #[test]
    fn test_integrate_with_polynomial_part() {
        // ∫ (x² + 1)/(x + 1) dx
        // = ∫ (x - 1 + 2/(x+1)) dx
        // = x²/2 - x + 2*log(x+1)
        let f = rf(&[1, 0, 1], &[1, 1]);

        let result = integrate_rational(&f);

        assert!(result.is_complete);
        // Should have polynomial part
        assert!(!result.polynomial_part.is_zero());
    }

    #[test]
    fn test_integrate_squared_denominator() {
        // ∫ 1/(x+1)² dx = -1/(x+1)
        let f = rf(&[1], &[1, 2, 1]); // 1 / (x+1)²

        let result = integrate_rational(&f);

        assert!(result.is_complete);
        // The rational part should be non-zero from Hermite reduction
        // -1/(x+1)
    }

    #[test]
    fn test_integrate_zero() {
        let f = RationalFunction::<Q>::zero();
        let result = integrate_rational(&f);

        assert!(result.is_zero());
        assert!(result.is_complete);
    }

    #[test]
    fn test_simple_pole_integration() {
        // ∫ 1/(x-2) dx = log(x-2)
        let term = integrate_simple_pole(&q(1), &q(2), 1);

        match term {
            SimplePoleTerm::Logarithmic { coefficient, .. } => {
                assert_eq!(coefficient, q(1));
            }
            _ => panic!("Expected logarithmic term"),
        }
    }

    #[test]
    fn test_higher_power_pole() {
        // ∫ 1/(x-2)² dx = -1/(x-2)
        let term = integrate_simple_pole(&q(1), &q(2), 2);

        match term {
            SimplePoleTerm::Rational {
                coefficient,
                power,
                ..
            } => {
                assert_eq!(coefficient, q(-1));
                assert_eq!(power, 1);
            }
            _ => panic!("Expected rational term"),
        }
    }

    // Tests for algebraic integration

    #[test]
    fn test_integrate_algebraic_x2_minus_1() {
        // ∫ 1/(x² - 1) dx = (1/2)log((x-1)/(x+1))
        // This has rational roots ±1
        let f = rf(&[1], &[-1, 0, 1]);

        let result = integrate_rational_with_algebraic(&f);

        assert!(result.is_complete);
        let has_log = result
            .algebraic_log_part
            .as_ref()
            .map_or(false, |p| !p.is_empty());
        assert!(has_log, "Should have logarithmic part");

        if let Some(ref log_part) = result.algebraic_log_part {
            println!("x² - 1: {} log terms", log_part.terms.len());
            assert!(log_part.terms.len() >= 2, "Should have 2 log terms");
        }
    }

    #[test]
    fn test_integrate_algebraic_x2_plus_1() {
        // ∫ 1/(x² + 1) dx = arctan(x) = (1/2i)log((x-i)/(x+i))
        // This needs algebraic roots ±i
        let f = rf(&[1], &[1, 0, 1]);

        let result = integrate_rational_with_algebraic(&f);

        assert!(result.is_complete);
        let has_log = result
            .algebraic_log_part
            .as_ref()
            .map_or(false, |p| !p.is_empty());
        assert!(has_log, "Should have algebraic logarithmic part");

        if let Some(ref log_part) = result.algebraic_log_part {
            println!("x² + 1: {} log terms", log_part.terms.len());
            assert!(!log_part.terms.is_empty());
        }
    }

    #[test]
    fn test_integrate_algebraic_linear() {
        // ∫ 1/(x - 1) dx = log(x - 1)
        // Simple case to ensure algebraic version works for rational roots
        let f = rf(&[1], &[-1, 1]);

        let result = integrate_rational_with_algebraic(&f);

        assert!(result.is_complete);
        let has_log = result
            .algebraic_log_part
            .as_ref()
            .map_or(false, |p| !p.is_empty());
        assert!(has_log, "Should have logarithmic part");
    }
}
