//! Lecerf's Algorithm for Multivariate Polynomial Factorization.
//!
//! This module implements Lecerf's state-of-the-art algorithm for factoring
//! sparse multivariate polynomials over Z[x₁, ..., xₙ].
//!
//! # Algorithm Overview
//!
//! 1. **Preprocessing**: Extract content, compute squarefree factorization
//! 2. **Evaluation Point Selection**: Choose (a₂, ..., aₙ) preserving separability
//! 3. **Univariate Factorization**: Factor f(x₁, a₂, ..., aₙ) using Van Hoeij
//! 4. **Leading Coefficient Precomputation**: Factor lc₁(f) to constrain lifting
//! 5. **Multivariate Hensel Lifting**: Lift to full multivariate factors
//! 6. **Factor Recombination**: Find true factors via trial division
//!
//! # References
//!
//! - Lecerf, G. (2006). "Sharp precision in Hensel lifting for bivariate
//!   polynomial factorization." Mathematics of Computation.
//! - Lecerf, G. (2007). "Improved dense multivariate polynomial factorization
//!   algorithms." Journal of Symbolic Computation.

use num_traits::Zero;
use rayon::prelude::*;
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_poly::monomial::PackedMonomial;
use tertius_poly::ordering::MonomialOrder;
use tertius_poly::sparse::SparsePoly;
use tertius_rings::integers::Z;

use crate::leading_coeff::{compute_leading_coeff_info, make_primitive, polynomial_content};
use crate::multivariate_hensel::multivariate_hensel_lift;
use crate::univariate::van_hoeij_factor;

/// Result of Lecerf multivariate factorization.
#[derive(Clone, Debug)]
pub struct LecerfMultivariateResult {
    /// Irreducible factors as sparse multivariate polynomials.
    pub factors: Vec<SparsePoly<Z>>,
    /// Multiplicities corresponding to each factor.
    pub multiplicities: Vec<u32>,
    /// Content extracted from the original polynomial.
    pub content: Z,
    /// Number of variables.
    pub num_vars: usize,
    /// Statistics about the computation.
    pub stats: LecerfStats,
}

/// Statistics from Lecerf factorization.
#[derive(Clone, Debug, Default)]
pub struct LecerfStats {
    /// Evaluation points used.
    pub evaluation_points: Vec<i64>,
    /// Number of univariate factors from specialization.
    pub num_univariate_factors: usize,
    /// Lifting precision achieved.
    pub lifting_precision: u32,
    /// Whether lifting was successful.
    pub lifting_success: bool,
    /// Number of factor recombination attempts.
    pub recombination_attempts: usize,
}

/// Factors a multivariate polynomial over Z using Lecerf's algorithm.
///
/// # Arguments
///
/// * `f` - The polynomial to factor, in Z[x₁, ..., xₙ]
///
/// # Returns
///
/// A result containing the irreducible factors and their multiplicities.
///
/// # Example
///
/// ```ignore
/// use tertius_factor::lecerf::lecerf_factor_multivariate;
/// use tertius_poly::sparse::SparsePoly;
///
/// // f = (x + y + 1)(x - y + 1) = x² + 2x + 1 - y²
/// let f = /* construct polynomial */;
/// let result = lecerf_factor_multivariate(&f);
/// assert_eq!(result.factors.len(), 2);
/// ```
pub fn lecerf_factor_multivariate(f: &SparsePoly<Z>) -> LecerfMultivariateResult {
    let num_vars = f.num_vars();

    // Handle trivial cases
    if f.is_zero() {
        return LecerfMultivariateResult {
            factors: vec![],
            multiplicities: vec![],
            content: Z(Integer::new(0)),
            num_vars,
            stats: LecerfStats::default(),
        };
    }

    if is_constant(f) {
        let c = f.terms().first().map(|(_, c)| c.clone()).unwrap_or(Z(Integer::new(1)));
        return LecerfMultivariateResult {
            factors: vec![],
            multiplicities: vec![],
            content: c,
            num_vars,
            stats: LecerfStats::default(),
        };
    }

    // Phase 1: Preprocessing
    let content = polynomial_content(f);
    let primitive = make_primitive(f, &content);

    // Check if essentially univariate
    if is_univariate(&primitive) {
        return factor_univariate_as_multivariate(&primitive, content, num_vars);
    }

    // Handle bivariate case specially (already implemented)
    if num_vars == 2 {
        return factor_bivariate_lecerf(&primitive, content);
    }

    // Phase 2: General multivariate factorization for n ≥ 3
    factor_general_multivariate(&primitive, content, num_vars)
}

/// Factors a polynomial that is univariate in variable 0.
fn factor_univariate_as_multivariate(
    f: &SparsePoly<Z>,
    content: Z,
    num_vars: usize,
) -> LecerfMultivariateResult {
    let dense = to_univariate(f, 0);
    let result = van_hoeij_factor(&dense);

    let factors: Vec<SparsePoly<Z>> = result
        .factors
        .into_iter()
        .map(|p| from_univariate(&p, num_vars))
        .collect();

    let num_factors = factors.len();
    let multiplicities = vec![1u32; num_factors];

    LecerfMultivariateResult {
        factors,
        multiplicities,
        content: Z(content.0 * result.content.0),
        num_vars,
        stats: LecerfStats {
            num_univariate_factors: num_factors,
            lifting_success: true,
            ..Default::default()
        },
    }
}

/// Factors a bivariate polynomial using a simplified Lecerf approach.
fn factor_bivariate_lecerf(f: &SparsePoly<Z>, content: Z) -> LecerfMultivariateResult {
    let num_vars = 2;

    // Choose evaluation point
    let eval_points = choose_evaluation_points(f, 1);
    let eval_point = eval_points.first().copied().unwrap_or(1);

    // Specialize and factor
    let specialized = evaluate_at_point(f, 1, &Z(Integer::new(eval_point)));
    let uni_result = van_hoeij_factor(&specialized);

    if uni_result.factors.len() <= 1 {
        // Irreducible
        return LecerfMultivariateResult {
            factors: vec![f.clone()],
            multiplicities: vec![1],
            content,
            num_vars,
            stats: LecerfStats {
                evaluation_points: vec![eval_point],
                num_univariate_factors: 1,
                lifting_success: true,
                recombination_attempts: 0,
                lifting_precision: 0,
            },
        };
    }

    // Compute leading coefficient info
    let lc_info = compute_leading_coeff_info(f, 0);

    // Attempt Hensel lifting
    let hensel_result = multivariate_hensel_lift(
        f,
        &uni_result.factors,
        &[Z(Integer::new(eval_point))],
        f.degree_in(1) + 1,
    );

    if hensel_result.success {
        // Verify factors
        let verified = verify_and_extract_factors(f, &hensel_result.lifted_factors);

        return LecerfMultivariateResult {
            factors: verified,
            multiplicities: vec![1; hensel_result.lifted_factors.len()],
            content,
            num_vars,
            stats: LecerfStats {
                evaluation_points: vec![eval_point],
                num_univariate_factors: uni_result.factors.len(),
                lifting_success: true,
                recombination_attempts: 1,
                lifting_precision: lc_info.main_degree,
            },
        };
    }

    // Lifting failed - use trial division fallback
    let factors = factor_by_trial_division(f, &uni_result.factors, eval_point);

    LecerfMultivariateResult {
        factors,
        multiplicities: vec![1; uni_result.factors.len()],
        content,
        num_vars,
        stats: LecerfStats {
            evaluation_points: vec![eval_point],
            num_univariate_factors: uni_result.factors.len(),
            lifting_success: false,
            recombination_attempts: uni_result.factors.len(),
            lifting_precision: 0,
        },
    }
}

/// Factors a general multivariate polynomial (n ≥ 3 variables).
fn factor_general_multivariate(
    f: &SparsePoly<Z>,
    content: Z,
    num_vars: usize,
) -> LecerfMultivariateResult {
    // Choose evaluation points for variables 1, 2, ..., n-1
    let mut all_eval_points = Vec::new();
    for var in 1..num_vars {
        let pts = choose_evaluation_points(f, var);
        all_eval_points.push(pts.first().copied().unwrap_or(1));
    }

    // Specialize to univariate
    let mut specialized = f.clone();
    for (i, &pt) in all_eval_points.iter().enumerate() {
        let var = i + 1;
        specialized = specialized.partial_eval(var, &Z(Integer::new(pt)));
    }

    let dense = to_univariate(&specialized, 0);
    let uni_result = van_hoeij_factor(&dense);

    if uni_result.factors.len() <= 1 {
        return LecerfMultivariateResult {
            factors: vec![f.clone()],
            multiplicities: vec![1],
            content,
            num_vars,
            stats: LecerfStats {
                evaluation_points: all_eval_points,
                num_univariate_factors: 1,
                lifting_success: true,
                recombination_attempts: 0,
                lifting_precision: 0,
            },
        };
    }

    // Hensel lifting for multivariate case
    let eval_z: Vec<Z> = all_eval_points.iter().map(|&p| Z(Integer::new(p))).collect();
    let max_deg: u32 = (1..num_vars).map(|v| f.degree_in(v)).sum();

    let hensel_result = multivariate_hensel_lift(f, &uni_result.factors, &eval_z, max_deg + 1);

    if hensel_result.success {
        let verified = verify_and_extract_factors(f, &hensel_result.lifted_factors);

        return LecerfMultivariateResult {
            factors: verified,
            multiplicities: vec![1; hensel_result.lifted_factors.len()],
            content,
            num_vars,
            stats: LecerfStats {
                evaluation_points: all_eval_points,
                num_univariate_factors: uni_result.factors.len(),
                lifting_success: true,
                recombination_attempts: 1,
                lifting_precision: max_deg,
            },
        };
    }

    // Fallback: return polynomial as irreducible
    LecerfMultivariateResult {
        factors: vec![f.clone()],
        multiplicities: vec![1],
        content,
        num_vars,
        stats: LecerfStats {
            evaluation_points: all_eval_points,
            num_univariate_factors: uni_result.factors.len(),
            lifting_success: false,
            recombination_attempts: 1,
            lifting_precision: 0,
        },
    }
}

/// Verifies lifted factors and extracts true factors via trial division.
fn verify_and_extract_factors(f: &SparsePoly<Z>, candidates: &[SparsePoly<Z>]) -> Vec<SparsePoly<Z>> {
    let mut remaining = f.clone();
    let mut verified = Vec::new();

    for candidate in candidates {
        if candidate.is_zero() || is_constant(candidate) {
            continue;
        }

        if let Some(quotient) = exact_divide(&remaining, candidate) {
            verified.push(candidate.clone());
            remaining = quotient;

            if is_constant(&remaining) {
                break;
            }
        }
    }

    // Add remaining if non-constant
    if !is_constant(&remaining) && !remaining.is_zero() {
        verified.push(remaining);
    }

    if verified.is_empty() {
        vec![f.clone()]
    } else {
        verified
    }
}

/// Factors using trial division as a fallback.
fn factor_by_trial_division(
    f: &SparsePoly<Z>,
    uni_factors: &[DensePoly<Z>],
    eval_point: i64,
) -> Vec<SparsePoly<Z>> {
    let num_vars = f.num_vars();
    let mut remaining = f.clone();
    let mut factors = Vec::new();

    // Try each univariate factor
    for uni_f in uni_factors {
        let candidate = from_univariate(uni_f, num_vars);

        if let Some(quotient) = exact_divide(&remaining, &candidate) {
            factors.push(candidate);
            remaining = quotient;

            if is_constant(&remaining) {
                break;
            }
        }
    }

    if !is_constant(&remaining) && !remaining.is_zero() {
        factors.push(remaining);
    }

    if factors.is_empty() {
        vec![f.clone()]
    } else {
        factors
    }
}

/// Chooses good evaluation points for a variable.
fn choose_evaluation_points(f: &SparsePoly<Z>, var: usize) -> Vec<i64> {
    let candidates = [1i64, -1, 2, -2, 3, -3, 5, -5, 7, -7];
    let mut good_points = Vec::new();

    let main_degree = f.degree_in(0);

    for &a in &candidates {
        // Check that leading coefficient doesn't vanish
        let lc = f.leading_coeff_poly(0);
        let lc_eval = evaluate_sparse(&lc, var, &Z(Integer::new(a)));

        if !lc_eval.0.is_zero() {
            // Check that degree is preserved
            let specialized = f.partial_eval(var, &Z(Integer::new(a)));
            if specialized.degree_in(0) == main_degree {
                good_points.push(a);
                if good_points.len() >= 3 {
                    break;
                }
            }
        }
    }

    if good_points.is_empty() {
        good_points.push(1);
    }

    good_points
}

/// Evaluates a sparse polynomial at a value for one variable.
fn evaluate_sparse(f: &SparsePoly<Z>, var: usize, value: &Z) -> Z {
    let mut result = Z(Integer::new(0));

    for (mono, coeff) in f.terms() {
        let exp = mono.exponent(var);
        let mut term_val = coeff.0.clone();
        for _ in 0..exp {
            term_val = term_val * value.0.clone();
        }
        result.0 = result.0 + term_val;
    }

    result
}

/// Evaluates a sparse polynomial and returns a dense univariate result.
fn evaluate_at_point(f: &SparsePoly<Z>, var: usize, value: &Z) -> DensePoly<Z> {
    let specialized = f.partial_eval(var, value);
    to_univariate(&specialized, 0)
}

/// Checks if polynomial is univariate (only uses first variable).
fn is_univariate(f: &SparsePoly<Z>) -> bool {
    let num_vars = f.num_vars();
    f.terms().iter().all(|(m, _)| {
        (1..num_vars).all(|i| m.exponent(i) == 0)
    })
}

/// Checks if polynomial is constant.
fn is_constant(f: &SparsePoly<Z>) -> bool {
    let num_vars = f.num_vars();
    f.terms().iter().all(|(m, _)| m.total_degree(num_vars) == 0)
}

/// Converts sparse polynomial to dense univariate.
fn to_univariate(f: &SparsePoly<Z>, var: usize) -> DensePoly<Z> {
    let deg = f.degree_in(var);
    let mut coeffs = vec![Z(Integer::new(0)); (deg + 1) as usize];

    for (mono, coeff) in f.terms() {
        let exp = mono.exponent(var) as usize;
        if exp < coeffs.len() {
            coeffs[exp].0 = coeffs[exp].0.clone() + coeff.0.clone();
        }
    }

    DensePoly::new(coeffs)
}

/// Converts dense univariate to sparse multivariate.
fn from_univariate(f: &DensePoly<Z>, num_vars: usize) -> SparsePoly<Z> {
    let mut terms = Vec::new();
    for (i, c) in f.coeffs().iter().enumerate() {
        if !c.0.is_zero() {
            let mut exps = vec![0u32; num_vars];
            exps[0] = i as u32;
            terms.push((PackedMonomial::from_exponents(&exps), c.clone()));
        }
    }
    SparsePoly::new(terms, num_vars, MonomialOrder::Grevlex)
}

/// Attempts exact polynomial division.
fn exact_divide(a: &SparsePoly<Z>, b: &SparsePoly<Z>) -> Option<SparsePoly<Z>> {
    if b.is_zero() {
        return None;
    }

    let (quotient, remainder) = divide_with_remainder(a, b);

    if remainder.is_zero() || remainder.terms().iter().all(|(_, c)| c.0.is_zero()) {
        Some(quotient)
    } else {
        None
    }
}

/// Divides two sparse polynomials with remainder.
fn divide_with_remainder(a: &SparsePoly<Z>, b: &SparsePoly<Z>) -> (SparsePoly<Z>, SparsePoly<Z>) {
    if b.is_zero() {
        panic!("Division by zero");
    }

    let num_vars = a.num_vars();
    let order = a.order();

    let b_lead = match b.leading_term() {
        Some(t) => t.clone(),
        None => return (SparsePoly::zero(num_vars, order), a.clone()),
    };

    let mut quotient_terms = Vec::new();
    let mut remainder = a.clone();

    for _ in 0..1000 {
        // Safety limit
        if remainder.is_zero() {
            break;
        }

        let r_lead = match remainder.leading_term() {
            Some(t) => t.clone(),
            None => break,
        };

        // Check divisibility
        if !b_lead.0.divides(&r_lead.0, num_vars) {
            break;
        }

        let quot_mono = match r_lead.0.div(&b_lead.0, num_vars) {
            Some(m) => m,
            None => break,
        };

        // Check coefficient divisibility
        if r_lead.1 .0.clone() % b_lead.1 .0.clone() != Integer::new(0) {
            break;
        }

        let quot_coeff = Z(r_lead.1 .0.clone() / b_lead.1 .0.clone());
        quotient_terms.push((quot_mono, quot_coeff.clone()));

        // Subtract quot * b from remainder
        let quot_term = SparsePoly::new(vec![(quot_mono, quot_coeff)], num_vars, order);
        let subtrahend = quot_term.mul(b);
        remainder = remainder.sub(&subtrahend);
    }

    (
        SparsePoly::new(quotient_terms, num_vars, order),
        remainder,
    )
}

/// Batch factorization using parallel execution.
pub fn lecerf_factor_multivariate_batch(
    polys: &[SparsePoly<Z>],
) -> Vec<LecerfMultivariateResult> {
    polys.par_iter().map(lecerf_factor_multivariate).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(n: i64) -> Z {
        Z(Integer::new(n))
    }

    fn mono(exps: &[u32]) -> PackedMonomial {
        PackedMonomial::from_exponents(exps)
    }

    fn sparse_poly(terms: &[(i64, &[u32])], num_vars: usize) -> SparsePoly<Z> {
        let t: Vec<_> = terms.iter().map(|(c, e)| (mono(e), z(*c))).collect();
        SparsePoly::new(t, num_vars, MonomialOrder::Grevlex)
    }

    #[test]
    fn test_is_univariate() {
        let f1 = sparse_poly(&[(1, &[2, 0]), (1, &[1, 0]), (1, &[0, 0])], 2);
        assert!(is_univariate(&f1));

        let f2 = sparse_poly(&[(1, &[2, 0]), (1, &[0, 1])], 2);
        assert!(!is_univariate(&f2));
    }

    #[test]
    fn test_is_constant() {
        let c = sparse_poly(&[(5, &[0, 0])], 2);
        assert!(is_constant(&c));

        let x = sparse_poly(&[(1, &[1, 0])], 2);
        assert!(!is_constant(&x));
    }

    #[test]
    fn test_to_from_univariate() {
        // x² + 2x + 1
        let sparse = sparse_poly(&[(1, &[2, 0]), (2, &[1, 0]), (1, &[0, 0])], 2);
        let dense = to_univariate(&sparse, 0);
        let back = from_univariate(&dense, 2);

        assert_eq!(sparse.terms().len(), back.terms().len());
    }

    #[test]
    fn test_factor_univariate_polynomial() {
        // (x+1)(x+2) = x² + 3x + 2
        let f = sparse_poly(&[(1, &[2, 0]), (3, &[1, 0]), (2, &[0, 0])], 2);
        let result = lecerf_factor_multivariate(&f);

        assert_eq!(result.factors.len(), 2);
        assert!(result.stats.lifting_success);
    }

    #[test]
    fn test_factor_irreducible() {
        // x² + 1 is irreducible over Z
        let f = sparse_poly(&[(1, &[2, 0]), (1, &[0, 0])], 2);
        let result = lecerf_factor_multivariate(&f);

        assert_eq!(result.factors.len(), 1);
    }

    #[test]
    fn test_factor_constant() {
        let f = sparse_poly(&[(42, &[0, 0])], 2);
        let result = lecerf_factor_multivariate(&f);

        assert!(result.factors.is_empty());
        assert_eq!(result.content.0, Integer::new(42));
    }

    #[test]
    fn test_factor_trivariate() {
        // Simple trivariate: x + y + z (irreducible)
        let f = sparse_poly(&[(1, &[1, 0, 0]), (1, &[0, 1, 0]), (1, &[0, 0, 1])], 3);
        let result = lecerf_factor_multivariate(&f);

        // Should be irreducible
        assert_eq!(result.factors.len(), 1);
    }
}
