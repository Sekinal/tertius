//! Leading coefficient utilities for multivariate factorization.
//!
//! Lecerf's algorithm requires careful handling of leading coefficients.
//! When factoring f(x₁,...,xₙ) viewed as f ∈ K[x₂,...,xₙ][x₁], the leading
//! coefficient lc₁(f) is a polynomial in the remaining variables.
//!
//! For each irreducible factor h of f, we have lc₁(h) | lc₁(f). This
//! constraint helps prune the search space during factor recombination.

use num_traits::{One, Zero};
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_poly::monomial::PackedMonomial;
use tertius_poly::ordering::MonomialOrder;
use tertius_poly::sparse::SparsePoly;
use tertius_rings::integers::Z;
use tertius_rings::traits::Ring;

use crate::univariate::van_hoeij_factor;

/// Information about the leading coefficient of a multivariate polynomial.
#[derive(Clone, Debug)]
pub struct LeadingCoeffInfo {
    /// The leading coefficient polynomial (in remaining variables).
    pub lc_polynomial: SparsePoly<Z>,
    /// Irreducible factors of the LC polynomial.
    pub lc_factors: Vec<SparsePoly<Z>>,
    /// All divisors that could appear as LC of a factor.
    pub divisor_candidates: Vec<SparsePoly<Z>>,
    /// The main variable used.
    pub main_var: usize,
    /// Degree in the main variable.
    pub main_degree: u32,
}

/// Computes leading coefficient information for multivariate factorization.
///
/// # Arguments
/// * `f` - The polynomial to analyze
/// * `main_var` - The main variable (typically 0 for x)
///
/// # Returns
/// Leading coefficient information including factors and divisor candidates.
pub fn compute_leading_coeff_info(f: &SparsePoly<Z>, main_var: usize) -> LeadingCoeffInfo {
    let main_degree = f.degree_in(main_var);
    let lc_poly = f.leading_coeff_poly(main_var);

    // If LC is constant, all factors have constant LC
    if is_constant(&lc_poly) {
        return LeadingCoeffInfo {
            lc_polynomial: lc_poly.clone(),
            lc_factors: vec![lc_poly.clone()],
            divisor_candidates: vec![lc_poly],
            main_var,
            main_degree,
        };
    }

    // Factor the LC polynomial recursively
    // For now, handle the univariate case when LC is in one variable
    let lc_factors = factor_lc_polynomial(&lc_poly);

    // Generate all divisor candidates (products of subsets of factors)
    let divisor_candidates = generate_divisors(&lc_poly, &lc_factors);

    LeadingCoeffInfo {
        lc_polynomial: lc_poly,
        lc_factors,
        divisor_candidates,
        main_var,
        main_degree,
    }
}

/// Attempts to factor the leading coefficient polynomial.
///
/// This is a simplified version that handles:
/// - Constant polynomials
/// - Univariate polynomials in any single variable
fn factor_lc_polynomial(lc: &SparsePoly<Z>) -> Vec<SparsePoly<Z>> {
    if is_constant(lc) {
        return vec![lc.clone()];
    }

    let num_vars = lc.num_vars();

    // Find which variable the LC is univariate in (if any)
    for var in 0..num_vars {
        if is_univariate_in_var(lc, var) {
            // Convert to dense and factor
            let dense = to_dense_in_var(lc, var);
            let result = van_hoeij_factor(&dense);

            // Convert factors back to sparse
            return result
                .factors
                .into_iter()
                .map(|f| from_dense_in_var(&f, var, num_vars, lc.order()))
                .collect();
        }
    }

    // If truly multivariate, return as single factor for now
    // Full implementation would recursively factor
    vec![lc.clone()]
}

/// Generates all divisors of a polynomial given its factorization.
///
/// For factorization p = c * f₁^e₁ * ... * fₖ^eₖ, generates all products
/// of subsets of factors that could be leading coefficients of factors.
fn generate_divisors(
    original: &SparsePoly<Z>,
    factors: &[SparsePoly<Z>],
) -> Vec<SparsePoly<Z>> {
    if factors.is_empty() {
        return vec![original.clone()];
    }

    let num_vars = original.num_vars();
    let order = original.order();

    // Start with the constant 1
    let one = SparsePoly::one(num_vars, order);
    let mut divisors = vec![one];

    // For each factor, extend divisors by multiplying
    for factor in factors {
        let mut new_divisors = divisors.clone();
        for existing in &divisors {
            let product = existing.mul(factor);
            // Check if product divides the original (approximately)
            if !new_divisors.iter().any(|d| sparse_eq(d, &product)) {
                new_divisors.push(product);
            }
        }
        divisors = new_divisors;
    }

    divisors
}

/// Checks compatibility between univariate factors and LC divisor candidates.
///
/// For each univariate factor g, finds which LC divisors could be the leading
/// coefficient of a multivariate lift of g.
///
/// # Arguments
/// * `lc_info` - Leading coefficient information
/// * `uni_factors` - Univariate factors from specialization
/// * `eval_point` - The evaluation point used for specialization
///
/// # Returns
/// For each univariate factor, a list of compatible LC candidates.
pub fn compatible_leading_coeffs(
    lc_info: &LeadingCoeffInfo,
    uni_factors: &[DensePoly<Z>],
    eval_point: &[Z],
) -> Vec<Vec<SparsePoly<Z>>> {
    let num_vars = lc_info.lc_polynomial.num_vars();

    uni_factors
        .iter()
        .map(|uni_f| {
            let uni_lc = if uni_f.is_zero() {
                Z(Integer::new(0))
            } else {
                uni_f.leading_coeff().clone()
            };

            // Find divisor candidates that evaluate to uni_lc at eval_point
            lc_info
                .divisor_candidates
                .iter()
                .filter(|d| {
                    let d_eval = evaluate_sparse(d, eval_point);
                    d_eval.0 == uni_lc.0 || is_unit_multiple(&d_eval, &uni_lc)
                })
                .cloned()
                .collect()
        })
        .collect()
}

/// Evaluates a sparse polynomial at a point (for secondary variables only).
fn evaluate_sparse(f: &SparsePoly<Z>, point: &[Z]) -> Z {
    let mut result = Z(Integer::new(0));

    for (mono, coeff) in f.terms() {
        let mut term_val = coeff.0.clone();
        for (i, val) in point.iter().enumerate() {
            let exp = mono.exponent(i);
            for _ in 0..exp {
                term_val = term_val * val.0.clone();
            }
        }
        result.0 = result.0 + term_val;
    }

    result
}

/// Checks if a divides b up to a unit (±1 in Z).
fn is_unit_multiple(a: &Z, b: &Z) -> bool {
    if a.0.is_zero() && b.0.is_zero() {
        return true;
    }
    if a.0.is_zero() || b.0.is_zero() {
        return false;
    }
    a.0.clone() == b.0.clone() || a.0.clone() == -b.0.clone()
}

/// Checks if a sparse polynomial is constant (degree 0 in all variables).
fn is_constant(f: &SparsePoly<Z>) -> bool {
    let num_vars = f.num_vars();
    f.terms().iter().all(|(m, _)| m.total_degree(num_vars) == 0)
}

/// Checks if polynomial is univariate in a specific variable.
fn is_univariate_in_var(f: &SparsePoly<Z>, var: usize) -> bool {
    let num_vars = f.num_vars();
    for (mono, _) in f.terms() {
        for i in 0..num_vars {
            if i != var && mono.exponent(i) > 0 {
                return false;
            }
        }
    }
    true
}

/// Converts sparse polynomial to dense in a specific variable.
fn to_dense_in_var(f: &SparsePoly<Z>, var: usize) -> DensePoly<Z> {
    let mut max_deg = 0u32;
    for (mono, _) in f.terms() {
        max_deg = max_deg.max(mono.exponent(var));
    }

    let mut coeffs = vec![Z(Integer::new(0)); (max_deg + 1) as usize];
    for (mono, coeff) in f.terms() {
        let exp = mono.exponent(var) as usize;
        coeffs[exp].0 = coeffs[exp].0.clone() + coeff.0.clone();
    }

    DensePoly::new(coeffs)
}

/// Converts dense polynomial back to sparse in a specific variable.
fn from_dense_in_var(
    f: &DensePoly<Z>,
    var: usize,
    num_vars: usize,
    order: MonomialOrder,
) -> SparsePoly<Z> {
    let mut terms = Vec::new();
    for (i, c) in f.coeffs().iter().enumerate() {
        if !c.0.is_zero() {
            let mut exps = vec![0u32; num_vars];
            exps[var] = i as u32;
            terms.push((PackedMonomial::from_exponents(&exps), c.clone()));
        }
    }
    SparsePoly::new(terms, num_vars, order)
}

/// Compares two sparse polynomials for equality.
fn sparse_eq(a: &SparsePoly<Z>, b: &SparsePoly<Z>) -> bool {
    if a.terms().len() != b.terms().len() {
        return false;
    }

    // Both are sorted, so we can compare term by term
    for ((m1, c1), (m2, c2)) in a.terms().iter().zip(b.terms().iter()) {
        if m1 != m2 || c1.0 != c2.0 {
            return false;
        }
    }
    true
}

/// Computes the content of a polynomial (GCD of coefficients).
pub fn polynomial_content(f: &SparsePoly<Z>) -> Z {
    f.terms()
        .iter()
        .fold(Z(Integer::new(0)), |acc, (_, c)| Z(gcd_int(&acc.0, &c.0)))
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

/// Makes a polynomial primitive by dividing out the content.
pub fn make_primitive(f: &SparsePoly<Z>, content: &Z) -> SparsePoly<Z> {
    if content.0.is_one() || content.0.is_zero() {
        return f.clone();
    }

    let terms: Vec<_> = f
        .terms()
        .iter()
        .map(|(m, c)| (*m, Z(c.0.clone() / content.0.clone())))
        .collect();

    SparsePoly::new(terms, f.num_vars(), f.order())
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
    fn test_is_constant() {
        let c = sparse_poly(&[(5, &[0, 0])], 2);
        assert!(is_constant(&c));

        let x = sparse_poly(&[(1, &[1, 0])], 2);
        assert!(!is_constant(&x));
    }

    #[test]
    fn test_is_univariate_in_var() {
        // x^2 + x + 1 is univariate in x (var 0)
        let f = sparse_poly(&[(1, &[2, 0]), (1, &[1, 0]), (1, &[0, 0])], 2);
        assert!(is_univariate_in_var(&f, 0));
        assert!(!is_univariate_in_var(&f, 1)); // not univariate in y

        // y^2 + y + 1 is univariate in y (var 1)
        let g = sparse_poly(&[(1, &[0, 2]), (1, &[0, 1]), (1, &[0, 0])], 2);
        assert!(!is_univariate_in_var(&g, 0));
        assert!(is_univariate_in_var(&g, 1));
    }

    #[test]
    fn test_to_dense_in_var() {
        // y^2 + 2y + 3
        let f = sparse_poly(&[(1, &[0, 2]), (2, &[0, 1]), (3, &[0, 0])], 2);
        let dense = to_dense_in_var(&f, 1);

        assert_eq!(dense.degree(), 2);
        assert_eq!(dense.coeffs()[0].0, Integer::new(3));
        assert_eq!(dense.coeffs()[1].0, Integer::new(2));
        assert_eq!(dense.coeffs()[2].0, Integer::new(1));
    }

    #[test]
    fn test_leading_coeff_info_constant() {
        // f = x^2 + x + 1 (LC w.r.t. x is 1)
        let f = sparse_poly(&[(1, &[2, 0]), (1, &[1, 0]), (1, &[0, 0])], 2);
        let info = compute_leading_coeff_info(&f, 0);

        assert!(is_constant(&info.lc_polynomial));
        assert_eq!(info.main_degree, 2);
    }

    #[test]
    fn test_leading_coeff_info_poly() {
        // f = (y+1)x^2 + x + 1, LC w.r.t. x is (y+1)
        let f = sparse_poly(&[(1, &[2, 1]), (1, &[2, 0]), (1, &[1, 0]), (1, &[0, 0])], 2);
        let info = compute_leading_coeff_info(&f, 0);

        // LC should be y + 1
        assert_eq!(info.lc_polynomial.terms().len(), 2);
        assert_eq!(info.main_degree, 2);
    }

    #[test]
    fn test_polynomial_content() {
        // f = 6x^2 + 4x + 2, content = 2
        let f = sparse_poly(&[(6, &[2, 0]), (4, &[1, 0]), (2, &[0, 0])], 2);
        let content = polynomial_content(&f);
        assert_eq!(content.0, Integer::new(2));
    }

    #[test]
    fn test_make_primitive() {
        // f = 6x^2 + 4x + 2 -> primitive = 3x^2 + 2x + 1
        let f = sparse_poly(&[(6, &[2, 0]), (4, &[1, 0]), (2, &[0, 0])], 2);
        let content = z(2);
        let prim = make_primitive(&f, &content);

        // Check the coefficients are divided by 2
        let terms: Vec<_> = prim.terms().iter().map(|(_, c)| c.0.clone()).collect();
        assert!(terms.contains(&Integer::new(3)));
        assert!(terms.contains(&Integer::new(2)));
        assert!(terms.contains(&Integer::new(1)));
    }

    #[test]
    fn test_evaluate_sparse() {
        // f = xy + y + 1, evaluate at y=2 (point = [_, 2])
        // But since we evaluate all vars, use point = [0, 2] for y only
        let f = sparse_poly(&[(1, &[1, 1]), (1, &[0, 1]), (1, &[0, 0])], 2);

        // f(0, 2) = 0*2 + 2 + 1 = 3
        let result = evaluate_sparse(&f, &[z(0), z(2)]);
        assert_eq!(result.0, Integer::new(3));
    }
}
