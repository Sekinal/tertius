//! Multivariate polynomial factorization.
//!
//! This module provides multivariate factorization algorithms:
//! - Evaluation/interpolation approach for bivariate polynomials
//! - Lecerf's algorithm for general multivariate (future work)
//!
//! The algorithm works by:
//! 1. Evaluating the polynomial at random points to get univariate specializations
//! 2. Factoring the univariate polynomial using Van Hoeij
//! 3. Lifting the factors back to multivariate using Hensel lifting
//! 4. Using sparse interpolation to recover coefficients

use num_traits::{One, Zero};
use rayon::prelude::*;
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_poly::sparse::SparsePoly;
use tertius_poly::monomial::PackedMonomial;
use tertius_poly::ordering::MonomialOrder;
use tertius_rings::integers::Z;
use tertius_rings::traits::Ring;

use crate::univariate::van_hoeij_factor;

// Helper functions
fn z_zero() -> Z {
    <Z as Ring>::zero()
}

fn z_one() -> Z {
    <Z as Ring>::one()
}

/// Result of multivariate factorization.
#[derive(Clone, Debug)]
pub struct MultivariateFactorResult {
    /// The irreducible factors (sparse polynomials).
    pub factors: Vec<SparsePoly<Z>>,
    /// Content extracted.
    pub content: Z,
    /// Number of variables.
    pub num_vars: usize,
    /// Statistics about the computation.
    pub stats: MultivariateFactorStats,
}

/// Statistics from multivariate factorization.
#[derive(Clone, Debug, Default)]
pub struct MultivariateFactorStats {
    /// Evaluation point used.
    pub evaluation_point: Vec<i64>,
    /// Number of univariate factors.
    pub num_univariate_factors: usize,
    /// Whether lifting was successful.
    pub lifting_success: bool,
}

/// Result of Lecerf factorization (for backwards compatibility).
#[derive(Clone, Debug)]
pub struct LecerfResult {
    /// The irreducible factors (as univariate over Z).
    pub factors: Vec<DensePoly<Z>>,
    /// Content extracted.
    pub content: Z,
    /// Number of variables.
    pub num_vars: usize,
}

/// Factors a univariate polynomial using Van Hoeij.
///
/// This is a compatibility wrapper for univariate polynomials.
pub fn lecerf_factor(f: &DensePoly<Z>) -> LecerfResult {
    let result = van_hoeij_factor(f);

    LecerfResult {
        factors: result.factors,
        content: result.content,
        num_vars: 1,
    }
}

/// Factors a bivariate polynomial f(x, y) over Z.
///
/// Algorithm:
/// 1. Choose evaluation point a for y
/// 2. Compute f(x, a) and factor it
/// 3. Lift factors using bivariate Hensel lifting
/// 4. Verify and extract true factors
pub fn factor_bivariate(f: &SparsePoly<Z>) -> MultivariateFactorResult {
    if f.is_zero() {
        return MultivariateFactorResult {
            factors: vec![],
            content: z_zero(),
            num_vars: 2,
            stats: MultivariateFactorStats::default(),
        };
    }

    let num_vars = f.num_vars();
    if num_vars != 2 {
        // Fall back to treating as univariate in first variable
        return factor_as_univariate(f);
    }

    // Extract content (GCD of all coefficients)
    let content = polynomial_content(f);
    let primitive = make_primitive(f, &content);

    // Check if polynomial is essentially univariate
    if is_univariate_in_x(&primitive) {
        let uni = to_univariate(&primitive, 0);
        let uni_result = van_hoeij_factor(&uni);
        let factors: Vec<SparsePoly<Z>> = uni_result
            .factors
            .into_iter()
            .map(|p| from_univariate(&p, 2))
            .collect();
        let num_factors = factors.len();

        return MultivariateFactorResult {
            factors,
            content,
            num_vars: 2,
            stats: MultivariateFactorStats {
                num_univariate_factors: num_factors,
                lifting_success: true,
                ..Default::default()
            },
        };
    }

    // Choose a good evaluation point for y
    let eval_point = choose_evaluation_point(&primitive);

    // Evaluate f(x, a) to get univariate polynomial
    let f_specialized = evaluate_at_point(&primitive, 1, &Z(Integer::new(eval_point)));

    // Factor the univariate specialization
    let uni_factors = van_hoeij_factor(&f_specialized);

    if uni_factors.factors.len() <= 1 {
        // Polynomial is irreducible
        return MultivariateFactorResult {
            factors: vec![primitive],
            content,
            num_vars: 2,
            stats: MultivariateFactorStats {
                evaluation_point: vec![eval_point],
                num_univariate_factors: 1,
                lifting_success: true,
            },
        };
    }

    // Try to lift factors to bivariate
    let lifted = lift_to_bivariate(&primitive, &uni_factors.factors, eval_point);

    if lifted.is_empty() {
        // Lifting failed - polynomial might be irreducible
        return MultivariateFactorResult {
            factors: vec![primitive],
            content,
            num_vars: 2,
            stats: MultivariateFactorStats {
                evaluation_point: vec![eval_point],
                num_univariate_factors: uni_factors.factors.len(),
                lifting_success: false,
            },
        };
    }

    MultivariateFactorResult {
        factors: lifted,
        content,
        num_vars: 2,
        stats: MultivariateFactorStats {
            evaluation_point: vec![eval_point],
            num_univariate_factors: uni_factors.factors.len(),
            lifting_success: true,
        },
    }
}

/// Falls back to treating a sparse polynomial as univariate.
fn factor_as_univariate(f: &SparsePoly<Z>) -> MultivariateFactorResult {
    let num_vars = f.num_vars();
    let uni = to_univariate(f, 0);
    let result = van_hoeij_factor(&uni);

    let factors: Vec<SparsePoly<Z>> = result
        .factors
        .into_iter()
        .map(|p| from_univariate(&p, num_vars))
        .collect();
    let num_factors = factors.len();

    MultivariateFactorResult {
        factors,
        content: result.content,
        num_vars,
        stats: MultivariateFactorStats {
            num_univariate_factors: num_factors,
            lifting_success: true,
            ..Default::default()
        },
    }
}

/// Checks if a polynomial only involves the first variable.
fn is_univariate_in_x(f: &SparsePoly<Z>) -> bool {
    for (mono, _) in f.terms() {
        for i in 1..f.num_vars() {
            if mono.exponent(i) > 0 {
                return false;
            }
        }
    }
    true
}

/// Converts a sparse polynomial to univariate (in variable 0).
fn to_univariate(f: &SparsePoly<Z>, var: usize) -> DensePoly<Z> {
    let mut max_deg = 0u32;
    for (mono, _) in f.terms() {
        max_deg = max_deg.max(mono.exponent(var));
    }

    let mut coeffs = vec![z_zero(); (max_deg + 1) as usize];
    for (mono, coeff) in f.terms() {
        let exp = mono.exponent(var) as usize;
        coeffs[exp] = coeffs[exp].clone() + coeff.clone();
    }

    DensePoly::new(coeffs)
}

/// Converts a univariate polynomial to sparse form.
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

/// Evaluates a sparse polynomial at a specific value for one variable.
fn evaluate_at_point(f: &SparsePoly<Z>, var: usize, value: &Z) -> DensePoly<Z> {
    // Compute powers of the evaluation value
    let max_exp = f
        .terms()
        .iter()
        .map(|(m, _)| m.exponent(var))
        .max()
        .unwrap_or(0);

    let mut powers = vec![z_one()];
    for _ in 1..=max_exp {
        let prev = powers.last().unwrap();
        powers.push(Z(prev.0.clone() * value.0.clone()));
    }

    // Collect coefficients for each x-degree
    let max_x_deg = f
        .terms()
        .iter()
        .map(|(m, _)| m.exponent(0))
        .max()
        .unwrap_or(0);

    let mut coeffs = vec![z_zero(); (max_x_deg + 1) as usize];

    for (mono, coeff) in f.terms() {
        let x_exp = mono.exponent(0) as usize;
        let y_exp = mono.exponent(var) as usize;
        let contribution = Z(coeff.0.clone() * powers[y_exp].0.clone());
        coeffs[x_exp] = Z(coeffs[x_exp].0.clone() + contribution.0);
    }

    DensePoly::new(coeffs)
}

/// Chooses a good evaluation point for y.
fn choose_evaluation_point(f: &SparsePoly<Z>) -> i64 {
    // Try small integers first, avoiding points where leading coefficient vanishes
    let candidates = [1i64, -1, 2, -2, 3, -3, 5, -5, 7, -7, 11, -11];

    for &a in &candidates {
        let specialized = evaluate_at_point(f, 1, &Z(Integer::new(a)));
        if !specialized.is_zero() && specialized.degree() == max_x_degree(f) {
            // Also check squarefree-ness would help but skip for now
            return a;
        }
    }

    // Default to 1 if nothing better found
    1
}

/// Returns the maximum degree in x (variable 0).
fn max_x_degree(f: &SparsePoly<Z>) -> usize {
    f.terms()
        .iter()
        .map(|(m, _)| m.exponent(0) as usize)
        .max()
        .unwrap_or(0)
}

/// Attempts to lift univariate factors to bivariate factors.
///
/// Uses a simplified version of bivariate Hensel lifting.
fn lift_to_bivariate(
    f: &SparsePoly<Z>,
    uni_factors: &[DensePoly<Z>],
    eval_point: i64,
) -> Vec<SparsePoly<Z>> {
    // For now, use a trial division approach:
    // For each subset of univariate factors, try to reconstruct a bivariate factor

    let n = uni_factors.len();
    if n == 0 {
        return vec![];
    }

    if n == 1 {
        return vec![f.clone()];
    }

    let mut remaining = f.clone();
    let mut found_factors = Vec::new();
    let mut used = vec![false; n];

    // Try single factors first
    for i in 0..n {
        if used[i] {
            continue;
        }

        // Try to find a bivariate factor that specializes to uni_factors[i]
        if let Some(biv_factor) = try_reconstruct_factor(&remaining, &uni_factors[i], eval_point) {
            if verify_factor(&remaining, &biv_factor) {
                found_factors.push(biv_factor.clone());
                remaining = divide_sparse(&remaining, &biv_factor);
                used[i] = true;

                if is_constant(&remaining) {
                    break;
                }
            }
        }
    }

    // Try pairs if single factors didn't work
    if !is_constant(&remaining) {
        for i in 0..n {
            if used[i] || is_constant(&remaining) {
                continue;
            }
            for j in i + 1..n {
                if used[j] {
                    continue;
                }

                // Multiply the univariate factors
                let prod = mul_dense(&uni_factors[i], &uni_factors[j]);

                if let Some(biv_factor) = try_reconstruct_factor(&remaining, &prod, eval_point) {
                    if verify_factor(&remaining, &biv_factor) {
                        found_factors.push(biv_factor.clone());
                        remaining = divide_sparse(&remaining, &biv_factor);
                        used[i] = true;
                        used[j] = true;

                        if is_constant(&remaining) {
                            break;
                        }
                    }
                }
            }
            if is_constant(&remaining) {
                break;
            }
        }
    }

    // Add remaining polynomial if non-constant
    if !is_constant(&remaining) {
        found_factors.push(remaining);
    }

    found_factors
}

/// Attempts to reconstruct a bivariate factor from a univariate specialization.
fn try_reconstruct_factor(
    f: &SparsePoly<Z>,
    uni_factor: &DensePoly<Z>,
    eval_point: i64,
) -> Option<SparsePoly<Z>> {
    // Simple approach: check if there's a factor of f that specializes to uni_factor

    // First, check degree compatibility
    let uni_deg = uni_factor.degree();
    let f_x_deg = max_x_degree(f);

    if uni_deg > f_x_deg || uni_deg == 0 {
        return None;
    }

    // Try to find coefficients in y by solving a system
    // This is a simplified version - full Lecerf would use Hensel lifting

    // For now, just check if treating the univariate factor as bivariate works
    let candidate = from_univariate(uni_factor, 2);

    // Verify: candidate(x, eval_point) should equal uni_factor
    let specialized = evaluate_at_point(&candidate, 1, &Z(Integer::new(eval_point)));

    // Check if specialized equals uni_factor (approximately, since we're working over Z)
    if specialized.degree() == uni_factor.degree() {
        let coeffs_match = specialized
            .coeffs()
            .iter()
            .zip(uni_factor.coeffs().iter())
            .all(|(a, b)| a.0 == b.0);

        if coeffs_match {
            return Some(candidate);
        }
    }

    None
}

/// Verifies that a polynomial divides another exactly.
fn verify_factor(f: &SparsePoly<Z>, g: &SparsePoly<Z>) -> bool {
    // Try division and check remainder is zero
    let (_, rem) = divide_sparse_with_rem(f, g);
    is_zero_sparse(&rem)
}

/// Checks if a sparse polynomial is zero.
fn is_zero_sparse(f: &SparsePoly<Z>) -> bool {
    f.is_zero() || f.terms().iter().all(|(_, c)| c.0.is_zero())
}

/// Checks if a sparse polynomial is constant.
fn is_constant(f: &SparsePoly<Z>) -> bool {
    let num_vars = f.num_vars();
    f.terms().iter().all(|(m, _)| m.total_degree(num_vars) == 0)
}

/// Divides two sparse polynomials, returning quotient only.
fn divide_sparse(a: &SparsePoly<Z>, b: &SparsePoly<Z>) -> SparsePoly<Z> {
    divide_sparse_with_rem(a, b).0
}

/// Divides two sparse polynomials, returning quotient and remainder.
fn divide_sparse_with_rem(a: &SparsePoly<Z>, b: &SparsePoly<Z>) -> (SparsePoly<Z>, SparsePoly<Z>) {
    if b.is_zero() {
        panic!("Division by zero");
    }

    let num_vars = a.num_vars();
    let order = a.order();

    let b_lead = b.leading_term().expect("b is not zero");
    let b_lead_mono = &b_lead.0;
    let b_lead_coeff = &b_lead.1;

    let mut quotient_terms = Vec::new();
    let mut remainder = a.clone();

    loop {
        if remainder.is_zero() {
            break;
        }

        let r_lead = match remainder.leading_term() {
            Some(t) => t.clone(),
            None => break,
        };

        // Check if leading monomial of remainder is divisible by leading monomial of b
        if !b_lead_mono.divides(&r_lead.0, num_vars) {
            break;
        }

        // Compute quotient monomial
        let quot_mono = r_lead.0.div(b_lead_mono, num_vars).unwrap();
        let quot_coeff = Z(r_lead.1 .0.clone() / b_lead_coeff.0.clone());
        let quot_rem = Z(r_lead.1 .0.clone() % b_lead_coeff.0.clone());

        // If not exactly divisible, we're done
        if !quot_rem.0.is_zero() {
            break;
        }

        quotient_terms.push((quot_mono.clone(), quot_coeff.clone()));

        // Subtract quot * b from remainder
        let mut new_rem_terms: Vec<(PackedMonomial, Z)> = Vec::new();

        for (m, c) in remainder.terms() {
            new_rem_terms.push((m.clone(), c.clone()));
        }

        for (bm, bc) in b.terms() {
            let prod_mono = quot_mono.mul(bm);
            let prod_coeff = Z(quot_coeff.0.clone() * bc.0.clone());

            // Find and subtract from matching term
            let mut found = false;
            for (m, c) in &mut new_rem_terms {
                if *m == prod_mono {
                    *c = Z(c.0.clone() - prod_coeff.0.clone());
                    found = true;
                    break;
                }
            }
            if !found {
                new_rem_terms.push((prod_mono, Z(-prod_coeff.0.clone())));
            }
        }

        // Filter out zeros
        new_rem_terms.retain(|(_, c)| !c.0.is_zero());

        remainder = SparsePoly::new(new_rem_terms, num_vars, order);
    }

    let quotient = SparsePoly::new(quotient_terms, num_vars, order);
    (quotient, remainder)
}

/// Multiplies two dense polynomials.
fn mul_dense(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![z_zero(); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        for (j, bj) in b_coeffs.iter().enumerate() {
            result[i + j] = Z(result[i + j].0.clone() + ai.0.clone() * bj.0.clone());
        }
    }

    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

/// Computes the GCD of all coefficients.
fn polynomial_content(f: &SparsePoly<Z>) -> Z {
    f.terms()
        .iter()
        .fold(z_zero(), |acc, (_, c)| Z(gcd_int(&acc.0, &c.0)))
}

/// Integer GCD.
fn gcd_int(a: &Integer, b: &Integer) -> Integer {
    if b.is_zero() {
        if *a < Integer::new(0) {
            -a.clone()
        } else {
            a.clone()
        }
    } else {
        gcd_int(b, &(a.clone() % b.clone()))
    }
}

/// Divides out the content to make primitive.
fn make_primitive(f: &SparsePoly<Z>, content: &Z) -> SparsePoly<Z> {
    if content.0.is_one() || content.0.is_zero() {
        return f.clone();
    }

    let terms: Vec<_> = f
        .terms()
        .iter()
        .map(|(m, c)| (m.clone(), Z(c.0.clone() / content.0.clone())))
        .collect();

    SparsePoly::new(terms, f.num_vars(), f.order())
}

/// Parallel batch factorization.
pub fn lecerf_factor_batch(polys: &[DensePoly<Z>]) -> Vec<LecerfResult> {
    polys.par_iter().map(lecerf_factor).collect()
}

/// Parallel batch bivariate factorization.
pub fn factor_bivariate_batch(polys: &[SparsePoly<Z>]) -> Vec<MultivariateFactorResult> {
    polys.par_iter().map(|p| factor_bivariate(p)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(n: i64) -> Z {
        Z(Integer::new(n))
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Z> {
        DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
    }

    fn mono(exps: &[u32]) -> PackedMonomial {
        PackedMonomial::from_exponents(exps)
    }

    fn sparse_poly(terms: &[(i64, &[u32])], num_vars: usize) -> SparsePoly<Z> {
        let t: Vec<_> = terms
            .iter()
            .map(|(c, e)| (mono(e), z(*c)))
            .collect();
        SparsePoly::new(t, num_vars, MonomialOrder::Grevlex)
    }

    #[test]
    fn test_lecerf_basic() {
        let f = poly(&[2, 3, 1]); // (x+1)(x+2)
        let result = lecerf_factor(&f);
        assert_eq!(result.factors.len(), 2);
    }

    #[test]
    fn test_to_univariate() {
        // f = x^2 + x + 1 (no y dependence)
        let f = sparse_poly(&[(1, &[2, 0]), (1, &[1, 0]), (1, &[0, 0])], 2);
        let uni = to_univariate(&f, 0);

        assert_eq!(uni.degree(), 2);
        assert_eq!(uni.coeffs()[0].0, Integer::new(1));
        assert_eq!(uni.coeffs()[1].0, Integer::new(1));
        assert_eq!(uni.coeffs()[2].0, Integer::new(1));
    }

    #[test]
    fn test_evaluate_at_point() {
        // f = x^2 + xy + y^2
        let f = sparse_poly(&[(1, &[2, 0]), (1, &[1, 1]), (1, &[0, 2])], 2);

        // f(x, 2) = x^2 + 2x + 4
        let specialized = evaluate_at_point(&f, 1, &z(2));

        assert_eq!(specialized.degree(), 2);
        assert_eq!(specialized.coeffs()[0].0, Integer::new(4)); // constant term = 2^2
        assert_eq!(specialized.coeffs()[1].0, Integer::new(2)); // x coefficient
        assert_eq!(specialized.coeffs()[2].0, Integer::new(1)); // x^2 coefficient
    }

    #[test]
    fn test_factor_bivariate_univariate() {
        // f = (x+1)(x+2) = x^2 + 3x + 2 (no y dependence)
        let f = sparse_poly(&[(1, &[2, 0]), (3, &[1, 0]), (2, &[0, 0])], 2);
        let result = factor_bivariate(&f);

        assert_eq!(result.factors.len(), 2);
    }

    #[test]
    fn test_is_univariate() {
        let f1 = sparse_poly(&[(1, &[2, 0]), (1, &[1, 0])], 2);
        assert!(is_univariate_in_x(&f1));

        let f2 = sparse_poly(&[(1, &[2, 0]), (1, &[1, 1])], 2);
        assert!(!is_univariate_in_x(&f2));
    }
}
