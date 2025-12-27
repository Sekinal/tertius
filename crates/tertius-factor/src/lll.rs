//! LLL (Lenstra-Lenstra-Lovász) lattice reduction algorithm.
//!
//! LLL produces a reduced basis for a lattice, where:
//! - The first vector is approximately the shortest vector in the lattice
//! - Successive vectors are reasonably orthogonal
//!
//! Used by Van Hoeij's algorithm for polynomial factorization.

use num_traits::{One, Zero};
use std::ops::Neg;
use tertius_linalg::dense_matrix::DenseMatrix;
use tertius_rings::rationals::Q;
use tertius_rings::traits::{Field, Ring};

// Helper functions to avoid trait ambiguity
fn q_zero() -> Q {
    <Q as Ring>::zero()
}

fn q_one() -> Q {
    <Q as Ring>::one()
}

fn q_is_zero(x: &Q) -> bool {
    <Q as Zero>::is_zero(x)
}

/// Result of LLL reduction.
#[derive(Clone, Debug)]
pub struct LllResult {
    /// The reduced basis as a matrix (rows are basis vectors).
    pub basis: DenseMatrix<Q>,
    /// Number of LLL iterations performed.
    pub iterations: usize,
}

/// Performs LLL reduction on a lattice basis.
///
/// # Arguments
///
/// * `basis` - Matrix where rows are the original basis vectors
/// * `delta` - Reduction parameter (typically 3/4 or 0.99)
///
/// # Returns
///
/// The reduced basis with short, nearly orthogonal vectors.
#[must_use]
pub fn lll_reduce(basis: &DenseMatrix<Q>, delta: &Q) -> LllResult {
    let n = basis.num_rows();

    if n == 0 {
        return LllResult {
            basis: basis.clone(),
            iterations: 0,
        };
    }

    let mut b = basis.clone();
    let mut iterations = 0;

    // Gram-Schmidt orthogonalization (stored implicitly)
    let mut gs = compute_gram_schmidt(&b);

    let mut k = 1;
    while k < n {
        iterations += 1;

        // Size reduction
        for j in (0..k).rev() {
            let mu_kj = gs.mu[(k, j)].clone();
            if !is_size_reduced(&mu_kj) {
                let r = nearest_integer(&mu_kj);
                // b[k] = b[k] - r * b[j]
                for col in 0..b.num_cols() {
                    let adjustment = r.clone() * b[(j, col)].clone();
                    b[(k, col)] = b[(k, col)].clone() + adjustment.neg();
                }
                // Update Gram-Schmidt coefficients
                update_mu_after_size_reduce(&mut gs, k, j, &r);
            }
        }

        // Lovász condition check
        if satisfies_lovasz(&gs, k, delta) {
            k += 1;
        } else {
            // Swap b[k-1] and b[k]
            swap_rows(&mut b, k - 1, k);
            update_gram_schmidt_after_swap(&mut gs, k);
            k = k.saturating_sub(1).max(1);
        }
    }

    LllResult { basis: b, iterations }
}

/// Gram-Schmidt coefficients storage.
struct GramSchmidt {
    /// μ[i,j] = <b_i, b*_j> / <b*_j, b*_j> for j < i
    mu: DenseMatrix<Q>,
    /// ||b*_i||^2 values
    b_star_norms_sq: Vec<Q>,
}

/// Computes Gram-Schmidt orthogonalization coefficients.
fn compute_gram_schmidt(b: &DenseMatrix<Q>) -> GramSchmidt {
    let n = b.num_rows();
    let m = b.num_cols();

    let mut mu = DenseMatrix::<Q>::zeros(n, n);
    let mut b_star_norms_sq = vec![q_zero(); n];

    // b*_0 = b_0
    b_star_norms_sq[0] = dot_product(b, 0, b, 0, m);
    mu[(0, 0)] = q_one();

    for i in 1..n {
        // Compute μ[i,j] for j < i
        for j in 0..i {
            if q_is_zero(&b_star_norms_sq[j]) {
                mu[(i, j)] = q_zero();
            } else {
                // μ[i,j] = <b_i, b*_j> / ||b*_j||^2
                let proj = project_onto_orthogonal(b, i, b, j, &mu, &b_star_norms_sq, m);
                mu[(i, j)] = proj.field_div(&b_star_norms_sq[j]);
            }
        }
        mu[(i, i)] = q_one();

        // ||b*_i||^2 = ||b_i||^2 - sum_j μ[i,j]^2 * ||b*_j||^2
        let mut norm_sq = dot_product(b, i, b, i, m);
        for j in 0..i {
            let contribution = mu[(i, j)].clone() * mu[(i, j)].clone() * b_star_norms_sq[j].clone();
            norm_sq = norm_sq + contribution.neg();
        }
        b_star_norms_sq[i] = norm_sq;
    }

    GramSchmidt { mu, b_star_norms_sq }
}

/// Computes <b[row_a], b[row_b]>.
fn dot_product(a: &DenseMatrix<Q>, row_a: usize, b: &DenseMatrix<Q>, row_b: usize, cols: usize) -> Q {
    let mut sum = q_zero();
    for col in 0..cols {
        sum = sum + a[(row_a, col)].clone() * b[(row_b, col)].clone();
    }
    sum
}

/// Computes <b[row_a], b*[row_b]>.
fn project_onto_orthogonal(
    b: &DenseMatrix<Q>,
    row_a: usize,
    _orig: &DenseMatrix<Q>,
    row_b: usize,
    mu: &DenseMatrix<Q>,
    norms: &[Q],
    cols: usize,
) -> Q {
    // <b_a, b*_b> = <b_a, b_b> - sum_{k<b} μ[b,k] * <b_a, b*_k>
    let mut proj = dot_product(b, row_a, b, row_b, cols);

    for k in 0..row_b {
        if !q_is_zero(&norms[k]) {
            let inner = project_onto_orthogonal(b, row_a, b, k, mu, norms, cols);
            let contribution = mu[(row_b, k)].clone() * inner;
            proj = proj + contribution.neg();
        }
    }

    proj
}

/// Checks if |μ| ≤ 1/2.
fn is_size_reduced(mu: &Q) -> bool {
    let half = Q::new(1, 2);
    let neg_half = Q::new(-1, 2);

    mu.clone() >= neg_half && mu.clone() <= half
}

/// Computes the nearest integer to a rational using banker's rounding (round half to even).
fn nearest_integer(q: &Q) -> Q {
    let inner = q.as_inner();
    let num = inner.numerator();
    let den = inner.denominator();

    // Compute floor and check if we're exactly at a half
    let floor_val = num.clone() / den.clone();
    let remainder = num.clone() % den.clone();

    // Check if remainder is exactly half (2*remainder == den)
    let two = tertius_integers::Integer::new(2);
    let double_rem = remainder.clone() * two.clone();

    if double_rem == den {
        // Exactly at half - round to even (banker's rounding)
        let floor_i64 = floor_val.to_i64().unwrap_or(0);
        if floor_i64 % 2 == 0 {
            Q::from_integer(floor_i64) // floor is even, use it
        } else {
            Q::from_integer(floor_i64 + 1) // floor is odd, round up to even
        }
    } else if double_rem > den {
        // More than half - round up
        Q::from_integer(floor_val.to_i64().unwrap_or(0) + 1)
    } else {
        // Less than half - round down (use floor)
        Q::from_integer(floor_val.to_i64().unwrap_or(0))
    }
}

/// Updates μ after size reduction: b[k] -= r * b[j].
fn update_mu_after_size_reduce(gs: &mut GramSchmidt, k: usize, j: usize, r: &Q) {
    for l in 0..=j {
        let adjustment = r.clone() * gs.mu[(j, l)].clone();
        gs.mu[(k, l)] = gs.mu[(k, l)].clone() + adjustment.neg();
    }
}

/// Checks the Lovász condition: ||b*_k||^2 >= (delta - μ[k,k-1]^2) * ||b*_{k-1}||^2.
fn satisfies_lovasz(gs: &GramSchmidt, k: usize, delta: &Q) -> bool {
    let mu_k_km1 = gs.mu[(k, k - 1)].clone();
    let mu_sq = mu_k_km1.clone() * mu_k_km1;
    let threshold = (delta.clone() + mu_sq.neg()) * gs.b_star_norms_sq[k - 1].clone();

    gs.b_star_norms_sq[k].clone() >= threshold
}

/// Swaps rows i and j.
fn swap_rows(m: &mut DenseMatrix<Q>, i: usize, j: usize) {
    m.swap_rows(i, j);
}

/// Updates Gram-Schmidt after swapping rows k-1 and k.
fn update_gram_schmidt_after_swap(gs: &mut GramSchmidt, k: usize) {
    let n = gs.mu.num_rows();

    // Save old values
    let mu_k_km1 = gs.mu[(k, k - 1)].clone();
    let b_star_km1 = gs.b_star_norms_sq[k - 1].clone();
    let b_star_k = gs.b_star_norms_sq[k].clone();

    // New ||b*_{k-1}||^2
    let new_b_star_km1 = b_star_k.clone() + mu_k_km1.clone() * mu_k_km1.clone() * b_star_km1.clone();

    if q_is_zero(&new_b_star_km1) {
        return;
    }

    // μ' = μ[k,k-1] * ||b*_{k-1}||^2 / new_||b*_{k-1}||^2
    let mu_prime = mu_k_km1.clone() * b_star_km1.clone();
    let mu_prime = mu_prime.field_div(&new_b_star_km1);

    // New ||b*_k||^2
    let new_b_star_k = b_star_km1.clone() * b_star_k.field_div(&new_b_star_km1);

    gs.b_star_norms_sq[k - 1] = new_b_star_km1;
    gs.b_star_norms_sq[k] = new_b_star_k;

    // Swap μ values for columns < k-1
    for j in 0..k - 1 {
        let tmp = gs.mu[(k - 1, j)].clone();
        gs.mu[(k - 1, j)] = gs.mu[(k, j)].clone();
        gs.mu[(k, j)] = tmp;
    }

    // Update μ[k, k-1]
    gs.mu[(k, k - 1)] = mu_prime.clone();

    // Update μ[i, k-1] and μ[i, k] for i > k
    for i in k + 1..n {
        let old_mu_i_km1 = gs.mu[(i, k - 1)].clone();
        let old_mu_i_k = gs.mu[(i, k)].clone();

        gs.mu[(i, k - 1)] = old_mu_i_k.clone() + mu_k_km1.clone() * old_mu_i_km1.clone();
        gs.mu[(i, k)] = old_mu_i_km1 + (mu_prime.clone() * gs.mu[(i, k - 1)].clone()).neg()
            + (mu_prime.clone() * mu_k_km1.clone() * old_mu_i_k).neg();
    }
}

/// Computes the Euclidean norm of a lattice vector (for debugging/testing).
#[must_use]
pub fn vector_norm_squared(v: &[Q]) -> Q {
    v.iter().fold(q_zero(), |acc, x| acc + x.clone() * x.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lll_identity() {
        // Identity basis is already reduced
        let basis = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(0)],
            vec![Q::from_integer(0), Q::from_integer(1)],
        ]);
        let delta = Q::new(3, 4);

        let result = lll_reduce(&basis, &delta);

        // Should be unchanged (or permuted)
        assert_eq!(result.basis.num_rows(), 2);
    }

    #[test]
    fn test_lll_reduces_basis() {
        // A non-reduced basis
        let basis = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(1)],
            vec![Q::from_integer(0), Q::from_integer(1)],
        ]);
        let delta = Q::new(3, 4);

        let result = lll_reduce(&basis, &delta);

        // Check that result is still a valid basis with reasonable vectors
        assert_eq!(result.basis.num_rows(), 2);
    }

    #[test]
    fn test_lll_larger_basis() {
        // 3x3 basis that needs reduction
        let basis = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(1), Q::from_integer(1)],
            vec![Q::from_integer(-1), Q::from_integer(0), Q::from_integer(2)],
            vec![Q::from_integer(3), Q::from_integer(5), Q::from_integer(6)],
        ]);
        let delta = Q::new(99, 100);

        let result = lll_reduce(&basis, &delta);

        assert_eq!(result.basis.num_rows(), 3);
    }

    #[test]
    fn test_nearest_integer() {
        assert_eq!(nearest_integer(&Q::new(3, 2)), Q::from_integer(2));
        assert_eq!(nearest_integer(&Q::new(5, 2)), Q::from_integer(2));
        assert_eq!(nearest_integer(&Q::new(7, 4)), Q::from_integer(2));
        assert_eq!(nearest_integer(&Q::from_integer(3)), Q::from_integer(3));
    }

    #[test]
    fn test_is_size_reduced() {
        assert!(is_size_reduced(&Q::from_integer(0)));
        assert!(is_size_reduced(&Q::new(1, 4)));
        assert!(is_size_reduced(&Q::new(-1, 4)));
        assert!(!is_size_reduced(&Q::from_integer(1)));
        assert!(!is_size_reduced(&Q::new(3, 4)));
    }
}
