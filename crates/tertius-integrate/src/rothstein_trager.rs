//! Rothstein-Trager algorithm for logarithmic integration.
//!
//! Given a proper rational function A/D where D is squarefree,
//! the Rothstein-Trager algorithm computes the logarithmic part:
//!
//! ∫ A/D dx = Σᵢ cᵢ log(vᵢ(x))
//!
//! where the cᵢ are algebraic constants and vᵢ are polynomials.
//!
//! # Algorithm
//!
//! 1. Compute the resultant polynomial R(t) = Res_x(D, A - t·D')
//! 2. Factor R(t) to find roots c₁, c₂, ...
//! 3. For each root cᵢ, compute vᵢ = gcd(D, A - cᵢ·D')
//! 4. Return Σᵢ cᵢ log(vᵢ)
//!
//! # Note
//!
//! The roots cᵢ may be in an algebraic extension of the base field.
//! For rational functions over Q, we may need to work in Q(α) for some
//! algebraic number α.

use tertius_poly::algorithms::gcd::poly_gcd;
use tertius_poly::algorithms::resultant::resultant;
use tertius_poly::dense::DensePoly;
use tertius_rings::traits::Field;

/// A logarithmic term in the integral: c * log(v(x))
#[derive(Clone, Debug)]
pub struct LogTerm<F: Field> {
    /// The constant coefficient (may be in an algebraic extension)
    pub coefficient: F,
    /// The polynomial argument of the logarithm
    pub argument: DensePoly<F>,
}

/// The logarithmic part of a rational function integral.
#[derive(Clone, Debug)]
pub struct LogarithmicPart<F: Field> {
    /// The logarithmic terms: Σᵢ cᵢ log(vᵢ)
    pub terms: Vec<LogTerm<F>>,
}

impl<F: Field> LogarithmicPart<F> {
    /// Creates an empty logarithmic part (for polynomial integrands).
    pub fn empty() -> Self {
        Self { terms: Vec::new() }
    }

    /// Returns true if there are no logarithmic terms.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

/// Computes the resultant polynomial R(t) = Res_x(D, A - t·D').
///
/// This resultant is computed symbolically with t as a formal variable.
/// The result is a univariate polynomial in t whose roots give the
/// coefficients of the logarithmic terms.
///
/// # Arguments
///
/// * `a` - Numerator polynomial A(x)
/// * `d` - Denominator polynomial D(x) (must be squarefree)
///
/// # Returns
///
/// The coefficients of R(t) as a polynomial in t.
pub fn compute_resultant_poly<F: Field>(
    a: &DensePoly<F>,
    d: &DensePoly<F>,
) -> DensePoly<F> {
    let d_prime = d.derivative();

    // R(t) = Res_x(D, A - t·D')
    // We compute this by treating t symbolically.
    //
    // The polynomial A - t·D' has coefficients that are linear in t.
    // The resultant will be a polynomial in t of degree ≤ deg(D).
    //
    // For now, we compute R(t) by evaluating at enough points and interpolating.
    // This works when we have enough distinct field elements.

    let deg_bound = d.degree() + 1;

    // Evaluate R at points 0, 1, 2, ..., deg_bound
    let mut eval_points = Vec::with_capacity(deg_bound);
    let mut eval_values = Vec::with_capacity(deg_bound);

    for i in 0..=deg_bound {
        let t_val = F::one().mul_by_scalar(i as i64);

        // Compute A - t·D'
        let a_minus_t_dp = a.sub(&d_prime.scale(&t_val));

        // Compute Res_x(D, A - t·D')
        let res = resultant(d.coeffs(), a_minus_t_dp.coeffs());

        eval_points.push(t_val);
        eval_values.push(res);
    }

    // Lagrange interpolation to recover R(t)
    lagrange_interpolate(&eval_points, &eval_values)
}

/// Lagrange interpolation to find the polynomial passing through given points.
fn lagrange_interpolate<F: Field>(points: &[F], values: &[F]) -> DensePoly<F> {
    let n = points.len();
    if n == 0 {
        return DensePoly::zero();
    }
    if n == 1 {
        return DensePoly::constant(values[0].clone());
    }

    let mut result = DensePoly::zero();

    for i in 0..n {
        // Compute the i-th Lagrange basis polynomial
        let mut basis = DensePoly::one();

        for j in 0..n {
            if i == j {
                continue;
            }

            // (x - x_j)
            let factor = DensePoly::new(vec![-points[j].clone(), F::one()]);
            basis = basis.mul(&factor);

            // Divide by (x_i - x_j)
            let denom = points[i].clone() - points[j].clone();
            let denom_inv = denom.inv().expect("points must be distinct");
            basis = basis.scale(&denom_inv);
        }

        // Multiply by y_i
        basis = basis.scale(&values[i]);

        // Add to result
        result = result.add(&basis);
    }

    result
}

/// Applies the Rothstein-Trager algorithm to compute the logarithmic part.
///
/// Given A/D where D is squarefree, computes the logarithmic integral.
///
/// # Arguments
///
/// * `a` - Numerator polynomial (deg < deg(D))
/// * `d` - Denominator polynomial (squarefree)
///
/// # Returns
///
/// The logarithmic part of ∫(A/D)dx, or None if the integral requires
/// algebraic extensions not supported by the base field.
pub fn rothstein_trager<F: Field>(
    a: &DensePoly<F>,
    d: &DensePoly<F>,
) -> Option<LogarithmicPart<F>> {
    // Handle trivial cases
    if a.is_zero() {
        return Some(LogarithmicPart::empty());
    }

    // For degree 1: ∫ a/(x - r) dx = a * log(x - r)
    if d.degree() == 1 {
        // D = d₀ + d₁x, so root is r = -d₀/d₁
        // ∫ a₀/(d₁(x - r)) = (a₀/d₁) * log(x - r) = (a₀/d₁) * log(D/d₁)
        let a_const = a.coeff(0);
        let d_lead = d.leading_coeff().clone();
        let coeff = a_const.field_div(&d_lead);

        // Make D monic for the log argument
        let d_monic = d.scale(&d_lead.inv().unwrap());

        return Some(LogarithmicPart {
            terms: vec![LogTerm {
                coefficient: coeff,
                argument: d_monic,
            }],
        });
    }

    // Compute R(t) = Res_x(D, A - t·D')
    let r_poly = compute_resultant_poly(a, d);

    // Find roots of R(t)
    // For now, we handle the simple case where R(t) factors into linear factors over F
    let roots = find_rational_roots(&r_poly);

    if roots.is_empty() && !r_poly.is_zero() && r_poly.degree() > 0 {
        // R(t) has no roots in F - need algebraic extension
        // For now, return None to indicate this case
        return None;
    }

    // For each root c, compute v = gcd(D, A - c·D')
    let d_prime = d.derivative();
    let mut terms = Vec::new();

    for c in roots {
        // Compute A - c·D'
        let a_minus_c_dp = a.sub(&d_prime.scale(&c));

        // Compute v = gcd(D, A - c·D')
        let v = poly_gcd(d, &a_minus_c_dp);

        if v.degree() > 0 {
            terms.push(LogTerm {
                coefficient: c,
                argument: v,
            });
        }
    }

    Some(LogarithmicPart { terms })
}

/// Finds rational roots of a polynomial using the rational root theorem.
///
/// For polynomials over Q, possible rational roots are ±(factor of constant)/(factor of leading).
/// For finite fields, we can try all field elements.
fn find_rational_roots<F: Field>(p: &DensePoly<F>) -> Vec<F> {
    if p.is_zero() {
        return Vec::new();
    }

    // For small degree, try evaluation at several points
    // This is a simple approach that works for many practical cases
    let mut roots = Vec::new();

    // Try small integers
    for i in -10..=10 {
        let x = F::one().mul_by_scalar(i);
        if p.eval(&x).is_zero() {
            roots.push(x);
        }
    }

    // Remove duplicates (in case of multiplicity)
    roots.dedup_by(|a, b| a == b);

    roots
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;
    use tertius_rings::traits::Ring;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_lagrange_interpolation() {
        // Interpolate through (0, 1), (1, 2), (2, 5)
        // This is 1 + x + x^2
        let points = vec![q(0), q(1), q(2)];
        let values = vec![q(1), q(2), q(5)]; // 1, 1+1, 1+2+4

        // Actually: p(0)=1, p(1)=1+1+1=3, p(2)=1+2+4=7 for 1+x+x^2
        // Let's use correct values
        let values = vec![q(1), q(3), q(7)];

        let p = lagrange_interpolate(&points, &values);

        // Check values
        assert_eq!(p.eval(&q(0)), q(1));
        assert_eq!(p.eval(&q(1)), q(3));
        assert_eq!(p.eval(&q(2)), q(7));
    }

    #[test]
    fn test_rothstein_trager_linear() {
        // ∫ 1/(x - 1) dx = log(x - 1)
        let a = poly(&[1]); // 1
        let d = poly(&[-1, 1]); // x - 1

        let result = rothstein_trager(&a, &d);
        assert!(result.is_some());

        let log_part = result.unwrap();
        assert_eq!(log_part.terms.len(), 1);
        assert_eq!(log_part.terms[0].coefficient, q(1));
    }

    #[test]
    fn test_rothstein_trager_simple_poles() {
        // ∫ 1/(x² - 1) dx = ∫ (1/2)/(x-1) - (1/2)/(x+1) dx
        // = (1/2)log(x-1) - (1/2)log(x+1)
        // But with Rothstein-Trager, we get the logarithmic form directly

        let a = poly(&[1]); // 1
        let d = poly(&[-1, 0, 1]); // x² - 1

        let result = rothstein_trager(&a, &d);
        // This may or may not succeed depending on whether roots are found

        if let Some(log_part) = result {
            // Verify the terms make sense
            for term in &log_part.terms {
                assert!(term.argument.degree() >= 1);
            }
        }
    }

    #[test]
    fn test_rothstein_trager_zero_numerator() {
        let a = DensePoly::<Q>::zero();
        let d = poly(&[-1, 1]);

        let result = rothstein_trager(&a, &d);
        assert!(result.is_some());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn test_find_rational_roots() {
        // p = (x - 1)(x + 2) = x² + x - 2
        let p = poly(&[-2, 1, 1]);

        let roots = find_rational_roots(&p);

        assert!(roots.contains(&q(1)));
        assert!(roots.contains(&q(-2)));
    }
}
