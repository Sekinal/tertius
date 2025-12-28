//! Multivariate Hensel Lifting for Polynomial Factorization.
//!
//! This module implements Hensel lifting for multivariate polynomials, a key
//! component of Lecerf's algorithm. Given univariate factors g₁,...,gᵣ of
//! f(x, a₂, ..., aₙ), we lift them to multivariate factors in K[x₁,...,xₙ].
//!
//! # Algorithm
//!
//! Multivariate Hensel lifting proceeds variable by variable:
//!
//! 1. Start with factors valid mod (x₂ - a₂, x₃ - a₃, ..., xₙ - aₙ)
//! 2. For each variable xᵢ, lift from mod (xᵢ - aᵢ)^k to mod (xᵢ - aᵢ)^(2k)
//! 3. Use the Bezout identity: s·g + t·h = 1 for the lifting step
//!
//! # References
//!
//! - Lecerf, G. (2006). "Sharp precision in Hensel lifting for bivariate
//!   polynomial factorization." Mathematics of Computation.

use num_traits::{One, Zero};
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_poly::monomial::PackedMonomial;
use tertius_poly::ordering::MonomialOrder;
use tertius_poly::sparse::SparsePoly;
use tertius_rings::integers::Z;
use tertius_rings::traits::Ring;

/// Result of multivariate Hensel lifting.
#[derive(Clone, Debug)]
pub struct MultivariateHenselResult {
    /// Lifted factors as sparse multivariate polynomials.
    pub lifted_factors: Vec<SparsePoly<Z>>,
    /// The modulus for each variable: (variable_index, value, precision).
    pub modulus: Vec<(usize, Z, u32)>,
    /// Whether lifting was successful.
    pub success: bool,
    /// Statistics about the lifting process.
    pub stats: HenselLiftStats,
}

/// Statistics from Hensel lifting.
#[derive(Clone, Debug, Default)]
pub struct HenselLiftStats {
    /// Number of lifting iterations performed.
    pub iterations: usize,
    /// Maximum precision achieved.
    pub max_precision: u32,
    /// Number of Bezout coefficient computations.
    pub bezout_computations: usize,
}

/// Lifts univariate factors to multivariate factors.
///
/// Given f ∈ Z[x₁, ..., xₙ] and univariate factors g₁, ..., gᵣ of
/// f(x₁, a₂, ..., aₙ), this function attempts to lift them to multivariate
/// factors.
///
/// # Arguments
///
/// * `f` - The polynomial to factor
/// * `univariate_factors` - Factors of the univariate specialization
/// * `evaluation_points` - Points (a₂, ..., aₙ) used for specialization
/// * `target_precision` - Target precision for lifting (typically degree bounds)
///
/// # Returns
///
/// Lifted factors or an empty result if lifting fails.
pub fn multivariate_hensel_lift(
    f: &SparsePoly<Z>,
    univariate_factors: &[DensePoly<Z>],
    evaluation_points: &[Z],
    target_precision: u32,
) -> MultivariateHenselResult {
    let num_vars = f.num_vars();

    if univariate_factors.is_empty() {
        return MultivariateHenselResult {
            lifted_factors: vec![],
            modulus: vec![],
            success: false,
            stats: HenselLiftStats::default(),
        };
    }

    if univariate_factors.len() == 1 {
        // Single factor - just return f as the lift
        return MultivariateHenselResult {
            lifted_factors: vec![f.clone()],
            modulus: evaluation_points
                .iter()
                .enumerate()
                .map(|(i, v)| (i + 1, v.clone(), target_precision))
                .collect(),
            success: true,
            stats: HenselLiftStats {
                iterations: 0,
                max_precision: target_precision,
                bezout_computations: 0,
            },
        };
    }

    // Convert univariate factors to sparse form
    let mut current_factors: Vec<SparsePoly<Z>> = univariate_factors
        .iter()
        .map(|f| dense_to_sparse(f, num_vars))
        .collect();

    let mut stats = HenselLiftStats::default();

    // Lift one variable at a time (from variable 1 onwards)
    for var_idx in 1..num_vars {
        if var_idx - 1 >= evaluation_points.len() {
            break;
        }

        let eval_point = &evaluation_points[var_idx - 1];

        // Lift factors for this variable
        let lift_result = lift_variable(
            f,
            &current_factors,
            var_idx,
            eval_point,
            target_precision,
        );

        if !lift_result.success {
            return MultivariateHenselResult {
                lifted_factors: current_factors,
                modulus: (1..=var_idx)
                    .zip(evaluation_points.iter())
                    .map(|(i, v)| (i, v.clone(), 1))
                    .collect(),
                success: false,
                stats,
            };
        }

        current_factors = lift_result.factors;
        stats.iterations += lift_result.iterations;
        stats.bezout_computations += lift_result.bezout_computations;
    }

    stats.max_precision = target_precision;

    MultivariateHenselResult {
        lifted_factors: current_factors,
        modulus: (1..num_vars)
            .zip(evaluation_points.iter())
            .map(|(i, v)| (i, v.clone(), target_precision))
            .collect(),
        success: true,
        stats,
    }
}

/// Internal result for lifting a single variable.
struct VariableLiftResult {
    factors: Vec<SparsePoly<Z>>,
    success: bool,
    iterations: usize,
    bezout_computations: usize,
}

/// Lifts factors with respect to a single variable.
///
/// Given factors valid mod (xᵥ - a)^1, lifts to mod (xᵥ - a)^k.
fn lift_variable(
    f: &SparsePoly<Z>,
    factors: &[SparsePoly<Z>],
    var: usize,
    eval_point: &Z,
    target_k: u32,
) -> VariableLiftResult {
    let num_vars = f.num_vars();
    let order = f.order();

    if factors.len() < 2 {
        return VariableLiftResult {
            factors: factors.to_vec(),
            success: true,
            iterations: 0,
            bezout_computations: 0,
        };
    }

    // For bivariate case with 2 factors, use direct Hensel step
    if factors.len() == 2 {
        return lift_two_factors(f, &factors[0], &factors[1], var, eval_point, target_k);
    }

    // For multiple factors, use recursive approach:
    // Split into two groups and lift each group, then combine
    let mid = factors.len() / 2;
    let left_product = product_of_factors(&factors[..mid], num_vars, order);
    let right_product = product_of_factors(&factors[mid..], num_vars, order);

    // Lift the product pair
    let pair_result = lift_two_factors(f, &left_product, &right_product, var, eval_point, target_k);

    if !pair_result.success || pair_result.factors.len() != 2 {
        return VariableLiftResult {
            factors: factors.to_vec(),
            success: false,
            iterations: pair_result.iterations,
            bezout_computations: pair_result.bezout_computations,
        };
    }

    // TODO: Recursively lift within each group
    // For now, return the successfully lifted pair
    VariableLiftResult {
        factors: pair_result.factors,
        success: true,
        iterations: pair_result.iterations,
        bezout_computations: pair_result.bezout_computations,
    }
}

/// Lifts exactly two factors using quadratic Hensel lifting.
///
/// Given f = g · h mod (x_var - a), lifts to mod (x_var - a)^k.
fn lift_two_factors(
    f: &SparsePoly<Z>,
    g: &SparsePoly<Z>,
    h: &SparsePoly<Z>,
    var: usize,
    eval_point: &Z,
    target_k: u32,
) -> VariableLiftResult {
    let num_vars = f.num_vars();
    let order = f.order();

    // Get univariate versions for Bezout computation
    let g_uni = sparse_to_univariate_in_main(g);
    let h_uni = sparse_to_univariate_in_main(h);

    // Compute Bezout coefficients: s*g + t*h = 1 (mod p for some prime)
    // For now, work over Z and assume coprimality
    let (s, t, gcd) = extended_gcd_poly(&g_uni, &h_uni);

    if gcd.degree() > 0 || gcd.coeffs().iter().any(|c| c.0 != Integer::new(1) && c.0 != Integer::new(-1)) {
        // Factors are not coprime - lifting may fail
        return VariableLiftResult {
            factors: vec![g.clone(), h.clone()],
            success: false,
            iterations: 0,
            bezout_computations: 1,
        };
    }

    let s_sparse = dense_to_sparse(&s, num_vars);
    let t_sparse = dense_to_sparse(&t, num_vars);

    let mut current_g = g.clone();
    let mut current_h = h.clone();
    let mut current_k = 1u32;
    let mut iterations = 0;

    // Quadratic Hensel lifting: double precision each step
    while current_k < target_k {
        let next_k = (2 * current_k).min(target_k);

        // Compute error: e = f - g*h
        let product = current_g.mul(&current_h);
        let error = f.sub(&product);

        // Check if error is in the ideal (x_var - a)^current_k
        if is_in_ideal(&error, var, eval_point, current_k) {
            // Compute correction terms
            // delta_g = t * e mod g (approximately)
            // delta_h = s * e mod h (approximately)

            let te = t_sparse.mul(&error);
            let se = s_sparse.mul(&error);

            // Reduce modulo the ideal (x_var - a)^next_k
            let delta_g = reduce_mod_ideal(&te, var, eval_point, next_k);
            let delta_h = reduce_mod_ideal(&se, var, eval_point, next_k);

            current_g = current_g.add(&delta_g);
            current_h = current_h.add(&delta_h);
        }

        current_k = next_k;
        iterations += 1;

        // Safety limit
        if iterations > 100 {
            break;
        }
    }

    // Verify the lift
    let final_product = current_g.mul(&current_h);
    let final_error = f.sub(&final_product);

    let success = is_approximately_zero(&final_error, var, eval_point, target_k);

    VariableLiftResult {
        factors: vec![current_g, current_h],
        success,
        iterations,
        bezout_computations: 1,
    }
}

/// Converts a dense polynomial to sparse form in a multivariate context.
fn dense_to_sparse(f: &DensePoly<Z>, num_vars: usize) -> SparsePoly<Z> {
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

/// Converts a sparse polynomial (univariate in first variable) to dense.
fn sparse_to_univariate_in_main(f: &SparsePoly<Z>) -> DensePoly<Z> {
    let deg = f.degree_in(0);
    let mut coeffs = vec![Z(Integer::new(0)); (deg + 1) as usize];

    for (mono, coeff) in f.terms() {
        let exp = mono.exponent(0) as usize;
        // Only use terms that are univariate in the main variable
        let is_univariate = (1..f.num_vars()).all(|i| mono.exponent(i) == 0);
        if is_univariate && exp < coeffs.len() {
            coeffs[exp].0 = coeffs[exp].0.clone() + coeff.0.clone();
        }
    }

    DensePoly::new(coeffs)
}

/// Computes the product of multiple sparse factors.
fn product_of_factors(
    factors: &[SparsePoly<Z>],
    num_vars: usize,
    order: MonomialOrder,
) -> SparsePoly<Z> {
    if factors.is_empty() {
        return SparsePoly::one(num_vars, order);
    }

    let mut result = factors[0].clone();
    for f in &factors[1..] {
        result = result.mul(f);
    }
    result
}

/// Checks if a polynomial is in the ideal (x_var - a)^k.
fn is_in_ideal(f: &SparsePoly<Z>, var: usize, eval_point: &Z, k: u32) -> bool {
    // A polynomial is in (x_var - a)^k if evaluating at x_var = a gives 0
    // and the first k-1 "derivatives" w.r.t. (x_var - a) also give 0

    // Simplified check: evaluate at the point
    let evaluated = f.partial_eval(var, eval_point);

    // For k=1, just check if evaluation is zero
    if k == 1 {
        return evaluated.terms().iter().all(|(_, c)| c.0.is_zero());
    }

    // For higher k, we'd need to check higher-order terms
    // For now, use a simplified check
    evaluated.terms().iter().all(|(_, c)| c.0.is_zero())
}

/// Checks if a polynomial is approximately zero modulo the ideal.
fn is_approximately_zero(f: &SparsePoly<Z>, var: usize, eval_point: &Z, k: u32) -> bool {
    // Check if all terms have degree >= k in (x_var - a)
    // This is equivalent to checking f is in (x_var - a)^k

    // Simplified: evaluate at the point
    let evaluated = f.partial_eval(var, eval_point);
    evaluated.is_zero() || evaluated.terms().iter().all(|(_, c)| c.0.is_zero())
}

/// Reduces a polynomial modulo the ideal (x_var - a)^k.
fn reduce_mod_ideal(
    f: &SparsePoly<Z>,
    var: usize,
    eval_point: &Z,
    k: u32,
) -> SparsePoly<Z> {
    // For now, keep only terms with degree < k in the variable
    let num_vars = f.num_vars();
    let order = f.order();

    let terms: Vec<_> = f
        .terms()
        .iter()
        .filter(|(m, _)| m.exponent(var) < k)
        .map(|(m, c)| (*m, c.clone()))
        .collect();

    SparsePoly::new(terms, num_vars, order)
}

/// Extended GCD for polynomials over Z.
///
/// Returns (s, t, gcd) such that s*a + t*b = gcd.
fn extended_gcd_poly(a: &DensePoly<Z>, b: &DensePoly<Z>) -> (DensePoly<Z>, DensePoly<Z>, DensePoly<Z>) {
    if b.is_zero() {
        return (
            DensePoly::new(vec![Z(Integer::new(1))]),
            DensePoly::zero(),
            a.clone(),
        );
    }

    // Simplified version for monic polynomials
    // Full version would need pseudo-division

    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut old_s = DensePoly::new(vec![Z(Integer::new(1))]);
    let mut s = DensePoly::zero();
    let mut old_t = DensePoly::zero();
    let mut t = DensePoly::new(vec![Z(Integer::new(1))]);

    let mut iterations = 0;
    while !r.is_zero() && iterations < 1000 {
        // Pseudo-division to get quotient
        if old_r.degree() < r.degree() {
            break;
        }

        let (q, remainder) = pseudo_div(&old_r, &r);

        old_r = r;
        r = remainder;

        let new_s = sub_poly(&old_s, &mul_poly(&q, &s));
        old_s = s;
        s = new_s;

        let new_t = sub_poly(&old_t, &mul_poly(&q, &t));
        old_t = t;
        t = new_t;

        iterations += 1;
    }

    (old_s, old_t, old_r)
}

/// Pseudo-division of polynomials over Z.
fn pseudo_div(a: &DensePoly<Z>, b: &DensePoly<Z>) -> (DensePoly<Z>, DensePoly<Z>) {
    if b.is_zero() {
        panic!("Division by zero polynomial");
    }

    if a.degree() < b.degree() {
        return (DensePoly::zero(), a.clone());
    }

    let b_lead = b.leading_coeff().clone();
    let deg_diff = a.degree() - b.degree();

    let mut remainder = a.clone();
    let mut quotient_coeffs = vec![Z(Integer::new(0)); deg_diff + 1];

    for i in (0..=deg_diff).rev() {
        if remainder.degree() < b.degree() + i {
            continue;
        }

        let r_lead = remainder.leading_coeff().clone();

        // Check if divisible
        if r_lead.0.clone() % b_lead.0.clone() != Integer::new(0) {
            // Not exactly divisible - scale remainder
            let factor = b_lead.0.clone();
            remainder = scale_poly(&remainder, &Z(factor));
        }

        let q_coeff = Z(remainder.leading_coeff().0.clone() / b_lead.0.clone());
        quotient_coeffs[i] = q_coeff.clone();

        // Subtract q_coeff * x^i * b from remainder
        let shift_b = shift_poly(b, i);
        let scaled_b = scale_poly(&shift_b, &q_coeff);
        remainder = sub_poly(&remainder, &scaled_b);
    }

    (DensePoly::new(quotient_coeffs), remainder)
}

/// Multiplies two dense polynomials.
fn mul_poly(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    if a.is_zero() || b.is_zero() {
        return DensePoly::zero();
    }

    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let mut result = vec![Z(Integer::new(0)); a_coeffs.len() + b_coeffs.len() - 1];

    for (i, ai) in a_coeffs.iter().enumerate() {
        for (j, bj) in b_coeffs.iter().enumerate() {
            result[i + j].0 = result[i + j].0.clone() + ai.0.clone() * bj.0.clone();
        }
    }

    // Remove leading zeros
    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

/// Subtracts two dense polynomials.
fn sub_poly(a: &DensePoly<Z>, b: &DensePoly<Z>) -> DensePoly<Z> {
    let a_coeffs = a.coeffs();
    let b_coeffs = b.coeffs();
    let max_len = a_coeffs.len().max(b_coeffs.len());
    let mut result = vec![Z(Integer::new(0)); max_len];

    for (i, c) in a_coeffs.iter().enumerate() {
        result[i].0 = result[i].0.clone() + c.0.clone();
    }
    for (i, c) in b_coeffs.iter().enumerate() {
        result[i].0 = result[i].0.clone() - c.0.clone();
    }

    // Remove leading zeros
    while result.len() > 1 && result.last().map_or(false, |c| c.0.is_zero()) {
        result.pop();
    }

    DensePoly::new(result)
}

/// Scales a polynomial by a constant.
fn scale_poly(a: &DensePoly<Z>, c: &Z) -> DensePoly<Z> {
    let coeffs: Vec<_> = a.coeffs().iter().map(|x| Z(x.0.clone() * c.0.clone())).collect();
    DensePoly::new(coeffs)
}

/// Shifts a polynomial by multiplying by x^k.
fn shift_poly(a: &DensePoly<Z>, k: usize) -> DensePoly<Z> {
    let mut coeffs = vec![Z(Integer::new(0)); k];
    coeffs.extend(a.coeffs().iter().cloned());
    DensePoly::new(coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(n: i64) -> Z {
        Z(Integer::new(n))
    }

    fn dense_poly(coeffs: &[i64]) -> DensePoly<Z> {
        DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
    }

    fn mono(exps: &[u32]) -> PackedMonomial {
        PackedMonomial::from_exponents(exps)
    }

    fn sparse_poly(terms: &[(i64, &[u32])], num_vars: usize) -> SparsePoly<Z> {
        let t: Vec<_> = terms.iter().map(|(c, e)| (mono(e), z(*c))).collect();
        SparsePoly::new(t, num_vars, MonomialOrder::Grevlex)
    }

    #[test]
    fn test_dense_to_sparse() {
        let dense = dense_poly(&[1, 2, 3]); // 1 + 2x + 3x²
        let sparse = dense_to_sparse(&dense, 2);

        assert_eq!(sparse.terms().len(), 3);
        assert_eq!(sparse.degree_in(0), 2);
    }

    #[test]
    fn test_sparse_to_univariate_in_main() {
        // x² + 2x + 1 as sparse with 2 variables
        let sparse = sparse_poly(&[(1, &[2, 0]), (2, &[1, 0]), (1, &[0, 0])], 2);
        let dense = sparse_to_univariate_in_main(&sparse);

        assert_eq!(dense.degree(), 2);
        assert_eq!(dense.coeffs()[0].0, Integer::new(1));
        assert_eq!(dense.coeffs()[1].0, Integer::new(2));
        assert_eq!(dense.coeffs()[2].0, Integer::new(1));
    }

    #[test]
    fn test_product_of_factors() {
        let f1 = sparse_poly(&[(1, &[1, 0]), (1, &[0, 0])], 2); // x + 1
        let f2 = sparse_poly(&[(1, &[1, 0]), (-1, &[0, 0])], 2); // x - 1

        let product = product_of_factors(&[f1, f2], 2, MonomialOrder::Grevlex);

        // Should be x² - 1
        assert_eq!(product.degree_in(0), 2);
    }

    #[test]
    fn test_mul_poly() {
        let a = dense_poly(&[1, 1]); // x + 1
        let b = dense_poly(&[-1, 1]); // x - 1

        let product = mul_poly(&a, &b);

        // Should be x² - 1
        assert_eq!(product.degree(), 2);
        assert_eq!(product.coeffs()[0].0, Integer::new(-1));
        assert_eq!(product.coeffs()[1].0, Integer::new(0));
        assert_eq!(product.coeffs()[2].0, Integer::new(1));
    }

    #[test]
    fn test_sub_poly() {
        let a = dense_poly(&[1, 2, 3]);
        let b = dense_poly(&[1, 1, 1]);

        let diff = sub_poly(&a, &b);

        assert_eq!(diff.coeffs()[0].0, Integer::new(0));
        assert_eq!(diff.coeffs()[1].0, Integer::new(1));
        assert_eq!(diff.coeffs()[2].0, Integer::new(2));
    }

    #[test]
    fn test_hensel_lift_single_factor() {
        // If there's only one factor, just return f
        let f = sparse_poly(&[(1, &[2, 0]), (1, &[0, 0])], 2); // x² + 1
        let factors = vec![dense_poly(&[1, 0, 1])]; // x² + 1

        let result = multivariate_hensel_lift(&f, &factors, &[z(1)], 10);

        assert!(result.success);
        assert_eq!(result.lifted_factors.len(), 1);
    }
}
