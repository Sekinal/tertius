//! Lecerf's algorithm for multivariate polynomial factorization.
//!
//! This module provides multivariate factorization via reduction to
//! univariate factorization. Full implementation deferred.

use rayon::prelude::*;
use tertius_poly::dense::DensePoly;
use tertius_rings::integers::Z;

use crate::univariate::van_hoeij_factor;

/// Result of multivariate factorization.
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
/// For true multivariate support, use the SparsePoly type directly.
/// This is a compatibility wrapper that works with DensePoly.
pub fn lecerf_factor(f: &DensePoly<Z>) -> LecerfResult {
    let result = van_hoeij_factor(f);

    LecerfResult {
        factors: result.factors,
        content: result.content,
        num_vars: 1,
    }
}

/// Parallel factorization.
pub fn lecerf_factor_batch(polys: &[DensePoly<Z>]) -> Vec<LecerfResult> {
    polys.par_iter().map(lecerf_factor).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_integers::Integer;

    fn z(n: i64) -> Z {
        Z(Integer::new(n))
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Z> {
        DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
    }

    #[test]
    fn test_lecerf_basic() {
        let f = poly(&[2, 3, 1]); // (x+1)(x+2)
        let result = lecerf_factor(&f);
        assert_eq!(result.factors.len(), 2);
    }
}
