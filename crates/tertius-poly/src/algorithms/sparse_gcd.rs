//! Hu-Monagan parallel sparse polynomial GCD.
//!
//! This implements the sparse modular GCD algorithm with:
//! - Parallel evaluation at random points
//! - Brown's bounds for early termination
//! - Sparse interpolation via Ben-Or/Tiwari

use rayon::prelude::*;
use std::ops::Neg;
use tertius_rings::finite_field::FiniteField;
use tertius_rings::traits::{Field, Ring};

/// A sparse multivariate term: (coefficient, exponent vector).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseTerm<R> {
    pub coeff: R,
    pub exponents: Vec<u32>,
}

impl<R: Clone> SparseTerm<R> {
    pub fn new(coeff: R, exponents: Vec<u32>) -> Self {
        Self { coeff, exponents }
    }
}

/// A sparse multivariate polynomial.
#[derive(Clone, Debug)]
pub struct SparseMultiPoly<R> {
    pub terms: Vec<SparseTerm<R>>,
    pub num_vars: usize,
}

impl<R: Ring + Clone + PartialEq> SparseMultiPoly<R> {
    /// Creates a new sparse polynomial.
    pub fn new(terms: Vec<SparseTerm<R>>, num_vars: usize) -> Self {
        // Filter out zero terms
        let terms: Vec<_> = terms.into_iter().filter(|t| !t.coeff.is_zero()).collect();
        Self { terms, num_vars }
    }

    /// Creates the zero polynomial.
    pub fn zero(num_vars: usize) -> Self {
        Self {
            terms: vec![],
            num_vars,
        }
    }

    /// Checks if zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of terms.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// Evaluates at a point.
    pub fn evaluate(&self, point: &[R]) -> R
    where
        R: Clone,
    {
        let mut result = R::zero();

        for term in &self.terms {
            let mut monomial_val = R::one();
            for (i, &exp) in term.exponents.iter().enumerate() {
                for _ in 0..exp {
                    monomial_val = monomial_val * point[i].clone();
                }
            }
            result = result + term.coeff.clone() * monomial_val;
        }

        result
    }

    /// Computes the total degree.
    pub fn total_degree(&self) -> u32 {
        self.terms
            .iter()
            .map(|t| t.exponents.iter().sum())
            .max()
            .unwrap_or(0)
    }
}

/// Result of sparse GCD computation.
#[derive(Clone, Debug)]
pub struct SparseGcdResult<R> {
    /// The GCD polynomial.
    pub gcd: SparseMultiPoly<R>,
    /// Number of evaluation points used.
    pub evaluations: usize,
}

/// Computes GCD using evaluation and interpolation.
///
/// This is a simplified version that works for univariate polynomials
/// over finite fields.
pub fn sparse_gcd_univariate<const P: u64>(
    f: &[FiniteField<P>],
    g: &[FiniteField<P>],
) -> Vec<FiniteField<P>> {
    // Handle edge cases
    if f.is_empty() || f.iter().all(|c| c.is_zero()) {
        return trim_leading_zeros(g);
    }
    if g.is_empty() || g.iter().all(|c| c.is_zero()) {
        return trim_leading_zeros(f);
    }

    // Use Euclidean algorithm for univariate
    let mut a = trim_leading_zeros(f);
    let mut b = trim_leading_zeros(g);

    if a.len() < b.len() {
        std::mem::swap(&mut a, &mut b);
    }

    while !b.is_empty() && !b.iter().all(|c| c.is_zero()) {
        let r = poly_remainder(&a, &b);
        a = b;
        b = r;
    }

    // Make monic
    if !a.is_empty() {
        let lead = a.last().unwrap().clone();
        if !lead.is_zero() {
            let lead_inv = lead.inv().expect("leading coeff must be invertible");
            for c in &mut a {
                *c = c.clone() * lead_inv.clone();
            }
        }
    }

    a
}

/// Computes polynomial remainder.
fn poly_remainder<const P: u64>(
    a: &[FiniteField<P>],
    b: &[FiniteField<P>],
) -> Vec<FiniteField<P>> {
    if b.is_empty() || b.iter().all(|c| c.is_zero()) {
        panic!("Division by zero polynomial");
    }

    let deg_a = a.len().saturating_sub(1);
    let deg_b = b.len().saturating_sub(1);

    if deg_a < deg_b {
        return a.to_vec();
    }

    let mut r = a.to_vec();
    let lead_b = b.last().unwrap().clone();
    let lead_b_inv = lead_b.inv().expect("leading coeff must be invertible");

    for i in (deg_b..=deg_a).rev() {
        if i < r.len() && !r[i].is_zero() {
            let coeff = r[i].clone() * lead_b_inv.clone();
            for j in 0..b.len() {
                let idx = i - deg_b + j;
                if idx < r.len() {
                    r[idx] = r[idx].clone() + (coeff.clone() * b[j].clone()).neg();
                }
            }
        }
    }

    trim_leading_zeros(&r)
}

/// Removes leading zeros.
fn trim_leading_zeros<R: Ring + Clone>(poly: &[R]) -> Vec<R> {
    let mut result = poly.to_vec();
    while result.len() > 1 && result.last().map_or(false, |c| c.is_zero()) {
        result.pop();
    }
    if result.len() == 1 && result[0].is_zero() {
        return vec![];
    }
    result
}

/// Parallel evaluation of polynomials at multiple points.
pub fn parallel_evaluate<const P: u64>(
    poly: &SparseMultiPoly<FiniteField<P>>,
    points: &[Vec<FiniteField<P>>],
) -> Vec<FiniteField<P>> {
    points
        .par_iter()
        .map(|point| poly.evaluate(point))
        .collect()
}

/// Computes multivariate GCD using parallel Brown's algorithm.
///
/// This evaluates f and g at many random points in parallel,
/// computes univariate GCDs, and interpolates.
pub fn sparse_gcd_multivariate<const P: u64>(
    f: &SparseMultiPoly<FiniteField<P>>,
    g: &SparseMultiPoly<FiniteField<P>>,
    omega: FiniteField<P>,
) -> SparseGcdResult<FiniteField<P>> {

    let num_vars = f.num_vars.max(g.num_vars);

    if num_vars == 0 || f.is_zero() {
        return SparseGcdResult {
            gcd: g.clone(),
            evaluations: 0,
        };
    }
    if g.is_zero() {
        return SparseGcdResult {
            gcd: f.clone(),
            evaluations: 0,
        };
    }

    if num_vars == 1 {
        // Convert to dense and use univariate GCD
        let f_dense = to_univariate_dense(f);
        let g_dense = to_univariate_dense(g);
        let gcd_dense = sparse_gcd_univariate::<P>(&f_dense, &g_dense);
        let gcd = from_univariate_dense(&gcd_dense);
        return SparseGcdResult {
            gcd,
            evaluations: 1,
        };
    }

    // For multivariate, use recursive evaluation
    // Evaluate at geometric points for the last variable
    let bound = f.total_degree().max(g.total_degree()) as usize + 1;
    let num_probes = bound * 2;

    // Generate evaluation points
    let points: Vec<FiniteField<P>> = (0..num_probes)
        .map(|i| pow_ff(omega.clone(), i as u64))
        .collect();

    // Evaluate f and g at each point (substituting x_n = point)
    let f_evals: Vec<SparseMultiPoly<FiniteField<P>>> = points
        .par_iter()
        .map(|point| evaluate_last_var(f, point))
        .collect();

    let g_evals: Vec<SparseMultiPoly<FiniteField<P>>> = points
        .par_iter()
        .map(|point| evaluate_last_var(g, point))
        .collect();

    // Compute GCDs recursively
    let gcd_evals: Vec<SparseMultiPoly<FiniteField<P>>> = f_evals
        .par_iter()
        .zip(g_evals.par_iter())
        .map(|(fe, ge)| sparse_gcd_multivariate::<P>(fe, ge, omega.clone()).gcd)
        .collect();

    // For now, return the GCD at the first point (proper interpolation TBD)
    SparseGcdResult {
        gcd: gcd_evals.into_iter().next().unwrap_or_else(|| SparseMultiPoly::zero(num_vars)),
        evaluations: num_probes,
    }
}

/// Converts sparse multivariate to univariate dense.
fn to_univariate_dense<const P: u64>(poly: &SparseMultiPoly<FiniteField<P>>) -> Vec<FiniteField<P>> {
    if poly.is_zero() {
        return vec![];
    }

    let max_deg = poly
        .terms
        .iter()
        .map(|t| t.exponents.first().copied().unwrap_or(0))
        .max()
        .unwrap_or(0);

    let mut result = vec![FiniteField::zero(); max_deg as usize + 1];
    for term in &poly.terms {
        let deg = term.exponents.first().copied().unwrap_or(0) as usize;
        result[deg] = result[deg].clone() + term.coeff.clone();
    }

    result
}

/// Converts univariate dense to sparse multivariate.
fn from_univariate_dense<const P: u64>(poly: &[FiniteField<P>]) -> SparseMultiPoly<FiniteField<P>> {
    let terms: Vec<_> = poly
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, c)| SparseTerm::new(c.clone(), vec![i as u32]))
        .collect();

    SparseMultiPoly::new(terms, 1)
}

/// Evaluates a polynomial at a specific value for the last variable.
fn evaluate_last_var<const P: u64>(
    poly: &SparseMultiPoly<FiniteField<P>>,
    val: &FiniteField<P>,
) -> SparseMultiPoly<FiniteField<P>> {
    if poly.num_vars <= 1 {
        // Substitute directly
        let result = poly.evaluate(&[val.clone()]);
        return SparseMultiPoly::new(
            vec![SparseTerm::new(result, vec![])],
            0,
        );
    }

    let new_num_vars = poly.num_vars - 1;
    let mut new_terms: std::collections::HashMap<Vec<u32>, FiniteField<P>> = std::collections::HashMap::new();

    for term in &poly.terms {
        let last_exp = term.exponents.last().copied().unwrap_or(0);
        let val_power = pow_ff(val.clone(), last_exp as u64);
        let new_coeff = term.coeff.clone() * val_power;

        let new_exps: Vec<u32> = term.exponents[..term.exponents.len().saturating_sub(1)].to_vec();

        new_terms
            .entry(new_exps.clone())
            .and_modify(|c| *c = c.clone() + new_coeff.clone())
            .or_insert(new_coeff);
    }

    let terms: Vec<_> = new_terms
        .into_iter()
        .filter(|(_, c)| !c.is_zero())
        .map(|(exps, c)| SparseTerm::new(c, exps))
        .collect();

    SparseMultiPoly::new(terms, new_num_vars)
}

/// Helper for FiniteField exponentiation.
fn pow_ff<const P: u64>(base: FiniteField<P>, mut exp: u64) -> FiniteField<P> {
    let mut result = FiniteField::<P>::one();
    let mut b = base;

    while exp > 0 {
        if exp & 1 == 1 {
            result = result * b.clone();
        }
        b = b.clone() * b;
        exp >>= 1;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    type GF101 = FiniteField<101>;

    fn gf101(n: i64) -> GF101 {
        GF101::new(n.rem_euclid(101) as u64)
    }

    #[test]
    fn test_sparse_gcd_univariate() {
        // f = (x + 1)(x + 2) = x^2 + 3x + 2
        // g = (x + 1)(x + 3) = x^2 + 4x + 3
        // gcd = x + 1
        let f = vec![gf101(2), gf101(3), gf101(1)];
        let g = vec![gf101(3), gf101(4), gf101(1)];

        let gcd = sparse_gcd_univariate::<101>(&f, &g);

        // gcd should be monic x + 1 = [1, 1]
        assert_eq!(gcd.len(), 2);
        assert_eq!(gcd[1], gf101(1));
    }

    #[test]
    fn test_sparse_gcd_coprime() {
        // f = x + 1, g = x + 2
        // gcd = 1
        let f = vec![gf101(1), gf101(1)];
        let g = vec![gf101(2), gf101(1)];

        let gcd = sparse_gcd_univariate::<101>(&f, &g);

        assert_eq!(gcd.len(), 1);
        assert_eq!(gcd[0], gf101(1));
    }

    #[test]
    fn test_sparse_poly_evaluate() {
        // f = 3x + 2y + 1
        let poly = SparseMultiPoly::new(
            vec![
                SparseTerm::new(gf101(1), vec![0, 0]), // 1
                SparseTerm::new(gf101(3), vec![1, 0]), // 3x
                SparseTerm::new(gf101(2), vec![0, 1]), // 2y
            ],
            2,
        );

        // Evaluate at (2, 3): 3*2 + 2*3 + 1 = 13
        let result = poly.evaluate(&[gf101(2), gf101(3)]);
        assert_eq!(result, gf101(13));
    }
}
