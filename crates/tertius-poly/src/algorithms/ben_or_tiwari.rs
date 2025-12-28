//! Ben-Or/Tiwari sparse polynomial interpolation.
//!
//! Given black-box access to a sparse polynomial f with t terms,
//! this algorithm recovers f using O(t) evaluations.
//!
//! Algorithm:
//! 1. Evaluate f at a geometric sequence: f(1), f(ω), f(ω²), ..., f(ω^{2t-1})
//! 2. Use Berlekamp-Massey to find the minimal polynomial whose roots are ω^{e_i}
//! 3. Find roots to recover the exponent vectors
//! 4. Solve a transposed Vandermonde system for coefficients

use crate::algorithms::berlekamp_massey::{berlekamp_massey, find_roots_exhaustive};
use rayon::prelude::*;
use std::ops::Neg;
use tertius_rings::finite_field::FiniteField;
use tertius_rings::traits::{Field, Ring};

/// Result of Ben-Or/Tiwari interpolation.
#[derive(Clone, Debug)]
pub struct BenOrTiwariResult<R> {
    /// Terms of the recovered polynomial: (coefficient, exponents).
    pub terms: Vec<(R, Vec<u32>)>,
    /// Number of evaluations used.
    pub num_evaluations: usize,
}

/// Interpolates a univariate sparse polynomial over a finite field.
///
/// # Arguments
/// - `evaluator`: Black-box evaluator f(x)
/// - `num_terms_bound`: Upper bound on the number of terms
/// - `omega`: Primitive element for geometric sequence
///
/// # Returns
/// The recovered polynomial as (coefficient, exponent) pairs.
pub fn ben_or_tiwari_univariate<const P: u64>(
    evaluator: impl Fn(FiniteField<P>) -> FiniteField<P> + Sync,
    num_terms_bound: usize,
    omega: FiniteField<P>,
) -> BenOrTiwariResult<FiniteField<P>> {
    let t = num_terms_bound;
    let num_evals = 2 * t;

    // Step 1: Evaluate at geometric sequence
    let evaluations: Vec<FiniteField<P>> = (0..num_evals)
        .into_par_iter()
        .map(|i| {
            let point = pow_ff(omega.clone(), i as u64);
            evaluator(point)
        })
        .collect();

    // Step 2: Berlekamp-Massey to find connection polynomial
    let bm_result = berlekamp_massey(&evaluations);
    let connection_poly = &bm_result.connection_poly;
    let actual_terms = bm_result.linear_complexity;

    if actual_terms == 0 {
        // Zero polynomial
        return BenOrTiwariResult {
            terms: vec![],
            num_evaluations: num_evals,
        };
    }

    // Step 3: Find roots of the connection polynomial
    // The roots are omega^{e_i} where e_i are the exponents
    let roots = find_roots_exhaustive(connection_poly);

    if roots.len() != actual_terms {
        // Failed to find all roots
        return BenOrTiwariResult {
            terms: vec![],
            num_evaluations: num_evals,
        };
    }

    // Step 4: Recover exponents from roots
    // omega^{e_i} = root => e_i = discrete_log(root, omega)
    let mut exponents = Vec::with_capacity(roots.len());
    for root in &roots {
        if let Some(exp) = discrete_log(&omega, root, P - 1) {
            exponents.push(exp as u32);
        }
    }

    if exponents.len() != actual_terms {
        return BenOrTiwariResult {
            terms: vec![],
            num_evaluations: num_evals,
        };
    }

    // Step 5: Solve transposed Vandermonde system for coefficients
    // V^T * c = y, where V_{ij} = (omega^{e_j})^i
    let coefficients = solve_transposed_vandermonde(&roots, &evaluations[..actual_terms]);

    // Build result
    let terms: Vec<_> = coefficients
        .into_iter()
        .zip(exponents.into_iter())
        .filter(|(c, _)| !c.is_zero())
        .map(|(c, e)| (c, vec![e]))
        .collect();

    BenOrTiwariResult {
        terms,
        num_evaluations: num_evals,
    }
}

/// Computes discrete logarithm: finds e such that base^e = target (mod p).
///
/// Uses baby-step giant-step for efficiency.
fn discrete_log<const P: u64>(
    base: &FiniteField<P>,
    target: &FiniteField<P>,
    order: u64,
) -> Option<u64> {
    if target.is_zero() {
        return None;
    }

    // Baby-step giant-step
    let m = ((order as f64).sqrt().ceil() as u64).max(1);

    // Baby steps: base^0, base^1, ..., base^{m-1}
    let mut baby_steps: std::collections::HashMap<u64, u64> = std::collections::HashMap::new();
    let mut power = FiniteField::<P>::one();
    for j in 0..m {
        baby_steps.insert(power.value(), j);
        power = power * base.clone();
    }

    // Giant step factor: base^{-m} = base^{order - m}
    let giant_factor = pow_ff(base.clone(), order - m % order);

    // Giant steps
    let mut gamma = target.clone();
    for i in 0..m {
        if let Some(&j) = baby_steps.get(&gamma.value()) {
            return Some(i * m + j);
        }
        gamma = gamma * giant_factor.clone();
    }

    None
}

/// Solves the transposed Vandermonde system V^T * c = y.
fn solve_transposed_vandermonde<const P: u64>(
    nodes: &[FiniteField<P>],
    values: &[FiniteField<P>],
) -> Vec<FiniteField<P>> {
    let n = nodes.len();
    if n == 0 {
        return vec![];
    }
    if n != values.len() {
        panic!("nodes and values must have same length");
    }

    // Use direct formula for small n, or iterative for larger
    if n == 1 {
        return vec![values[0].clone()];
    }

    // For the transposed Vandermonde, we need to solve:
    // sum_j c_j * nodes[j]^i = values[i]
    //
    // This can be done via Lagrange interpolation in a different form,
    // or by solving the linear system directly.

    // Simple Gaussian elimination for now (O(n^3))
    let mut matrix: Vec<Vec<FiniteField<P>>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            // V^T_{ij} = nodes[j]^i
            row.push(pow_ff(nodes[j].clone(), i as u64));
        }
        matrix.push(row);
    }

    // Augmented matrix
    let mut aug: Vec<Vec<FiniteField<P>>> = matrix
        .into_iter()
        .zip(values.iter())
        .map(|(row, val)| {
            let mut r = row;
            r.push(val.clone());
            r
        })
        .collect();

    // Forward elimination
    for col in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for row in col..n {
            if !aug[row][col].is_zero() {
                pivot_row = Some(row);
                break;
            }
        }

        let pivot_row = match pivot_row {
            Some(r) => r,
            None => continue,
        };

        // Swap rows
        aug.swap(col, pivot_row);

        // Scale pivot row
        let pivot = aug[col][col].clone();
        let pivot_inv = pivot.inv().expect("pivot must be invertible");
        for j in col..=n {
            aug[col][j] = aug[col][j].clone() * pivot_inv.clone();
        }

        // Eliminate
        for row in 0..n {
            if row != col {
                let factor = aug[row][col].clone();
                for j in col..=n {
                    let adj = factor.clone() * aug[col][j].clone();
                    aug[row][j] = aug[row][j].clone() + adj.neg();
                }
            }
        }
    }

    // Extract solution
    aug.into_iter().map(|row| row[n].clone()).collect()
}

/// Helper function for exponentiation.
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

    fn gf101(n: u64) -> GF101 {
        GF101::new(n % 101)
    }

    #[test]
    fn test_discrete_log() {
        // 2 is a primitive root mod 101
        let base = gf101(2);

        // 2^10 mod 101 = 1024 mod 101 = 14
        let target = gf101(14);
        let log = discrete_log(&base, &target, 100);
        assert_eq!(log, Some(10));
    }

    #[test]
    fn test_transposed_vandermonde() {
        // Simple 2x2 case
        let nodes = vec![gf101(2), gf101(3)];
        // V^T = [[1, 1], [2, 3]]
        // Solve [1, 1; 2, 3] * c = [5, 13]
        // c = [2, 3]
        let values = vec![gf101(5), gf101(13)];
        let coeffs = solve_transposed_vandermonde(&nodes, &values);

        assert_eq!(coeffs.len(), 2);
        assert_eq!(coeffs[0], gf101(2));
        assert_eq!(coeffs[1], gf101(3));
    }
}
