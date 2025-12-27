//! Parallel linear algebra operations.
//!
//! This module provides parallelized versions of common linear algebra
//! operations using rayon for work-stealing parallelism.

use num_traits::One;
use rayon::prelude::*;

use tertius_rings::traits::{Field, Ring};

use crate::dense_matrix::DenseMatrix;
use crate::sparse_matrix::CsrMatrix;

/// Configuration for parallel Gaussian elimination.
#[derive(Clone, Debug)]
pub struct ParallelConfig {
    /// Minimum matrix size to enable parallelism.
    pub parallel_threshold: usize,
    /// Block size for tiled algorithms.
    pub block_size: usize,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        Self {
            parallel_threshold: 64,
            block_size: 32,
        }
    }
}

/// Parallel row reduction for dense matrices over a field.
///
/// Uses tiled Gaussian elimination with rayon for parallelism.
pub fn parallel_row_reduce<R: Field + Clone + Send + Sync + One>(
    matrix: &mut DenseMatrix<R>,
    config: &ParallelConfig,
) -> usize {
    let num_rows = matrix.num_rows();
    let num_cols = matrix.num_cols();

    if num_rows < config.parallel_threshold {
        // Fall back to sequential for small matrices
        let (_, _, rank) = matrix.row_echelon();
        return rank;
    }

    let mut pivot_row = 0;
    let mut pivot_col = 0;

    while pivot_row < num_rows && pivot_col < num_cols {
        // Find pivot (sequential - typically fast)
        let mut max_row = None;
        for row in pivot_row..num_rows {
            if !matrix[(row, pivot_col)].is_zero() {
                max_row = Some(row);
                break;
            }
        }

        let Some(max_row) = max_row else {
            pivot_col += 1;
            continue;
        };

        // Swap rows if needed
        if max_row != pivot_row {
            matrix.swap_rows(pivot_row, max_row);
        }

        // Scale pivot row
        let pivot_val = matrix[(pivot_row, pivot_col)].clone();
        if let Some(inv) = pivot_val.inv() {
            matrix.scale_row(pivot_row, &inv);
        }

        // Parallel elimination of rows below pivot
        let pivot_row_data: Vec<R> = matrix.row(pivot_row).to_vec();

        // Collect rows to update
        let rows_to_update: Vec<usize> = (pivot_row + 1..num_rows)
            .filter(|&row| !matrix[(row, pivot_col)].is_zero())
            .collect();

        // Parallel update using work-stealing
        let updates: Vec<(usize, Vec<R>)> = rows_to_update
            .par_iter()
            .map(|&row| {
                let factor = matrix[(row, pivot_col)].clone().neg();
                let new_row: Vec<R> = matrix
                    .row(row)
                    .iter()
                    .zip(pivot_row_data.iter())
                    .map(|(a, b)| a.clone() + b.clone() * factor.clone())
                    .collect();
                (row, new_row)
            })
            .collect();

        // Apply updates (sequential - modifying shared state)
        for (row, new_row) in updates {
            for (col, val) in new_row.into_iter().enumerate() {
                matrix[(row, col)] = val;
            }
        }

        pivot_row += 1;
        pivot_col += 1;
    }

    pivot_row // rank
}

/// Parallel back-substitution for RREF.
pub fn parallel_back_substitute<R: Field + Clone + Send + Sync + One>(
    matrix: &mut DenseMatrix<R>,
    rank: usize,
) {
    let num_cols = matrix.num_cols();

    for pivot_row in (0..rank).rev() {
        // Find pivot column
        let pivot_col = (0..num_cols)
            .find(|&col| !matrix[(pivot_row, col)].is_zero())
            .unwrap_or(num_cols);

        if pivot_col >= num_cols {
            continue;
        }

        let pivot_row_data: Vec<R> = matrix.row(pivot_row).to_vec();

        // Rows to update (above pivot)
        let rows_to_update: Vec<usize> = (0..pivot_row)
            .filter(|&row| !matrix[(row, pivot_col)].is_zero())
            .collect();

        // Parallel update
        let updates: Vec<(usize, Vec<R>)> = rows_to_update
            .par_iter()
            .map(|&row| {
                let factor = matrix[(row, pivot_col)].clone().neg();
                let new_row: Vec<R> = matrix
                    .row(row)
                    .iter()
                    .zip(pivot_row_data.iter())
                    .map(|(a, b)| a.clone() + b.clone() * factor.clone())
                    .collect();
                (row, new_row)
            })
            .collect();

        for (row, new_row) in updates {
            for (col, val) in new_row.into_iter().enumerate() {
                matrix[(row, col)] = val;
            }
        }
    }
}

/// Parallel RREF computation.
pub fn parallel_rref<R: Field + Clone + Send + Sync + One>(
    matrix: &mut DenseMatrix<R>,
    config: &ParallelConfig,
) -> usize {
    let rank = parallel_row_reduce(matrix, config);
    parallel_back_substitute(matrix, rank);
    rank
}

/// Parallel sparse matrix-vector multiplication.
///
/// This is a simple wrapper around CsrMatrix::spmv_parallel.
pub fn parallel_spmv<R: Ring + Clone + Send + Sync>(matrix: &CsrMatrix<R>, x: &[R]) -> Vec<R> {
    matrix.spmv_parallel(x)
}

/// Parallel dot product of two vectors.
pub fn parallel_dot<R: Ring + Clone + Send + Sync>(a: &[R], b: &[R]) -> R {
    assert_eq!(a.len(), b.len());

    a.par_iter()
        .zip(b.par_iter())
        .map(|(ai, bi)| ai.clone() * bi.clone())
        .reduce(R::zero, |acc, x| acc + x)
}

/// Parallel vector addition: c = a + b.
pub fn parallel_add<R: Ring + Clone + Send + Sync>(a: &[R], b: &[R]) -> Vec<R> {
    assert_eq!(a.len(), b.len());

    a.par_iter()
        .zip(b.par_iter())
        .map(|(ai, bi)| ai.clone() + bi.clone())
        .collect()
}

/// Parallel scalar-vector multiplication: c = alpha * a.
pub fn parallel_scale<R: Ring + Clone + Send + Sync>(alpha: &R, a: &[R]) -> Vec<R> {
    a.par_iter().map(|ai| alpha.clone() * ai.clone()).collect()
}

/// Parallel axpy: y = alpha * x + y.
pub fn parallel_axpy<R: Ring + Clone + Send + Sync>(alpha: &R, x: &[R], y: &mut [R]) {
    assert_eq!(x.len(), y.len());

    y.par_iter_mut().zip(x.par_iter()).for_each(|(yi, xi)| {
        *yi = yi.clone() + alpha.clone() * xi.clone();
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_parallel_row_reduce() {
        let mut m = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)],
            vec![Q::from_integer(4), Q::from_integer(5), Q::from_integer(6)],
            vec![Q::from_integer(7), Q::from_integer(8), Q::from_integer(9)],
        ]);

        let config = ParallelConfig {
            parallel_threshold: 1, // Force parallel path for testing
            block_size: 2,
        };

        let rank = parallel_row_reduce(&mut m, &config);
        // Matrix has rank 2 (rows are linearly dependent)
        assert_eq!(rank, 2);
    }

    #[test]
    fn test_parallel_dot() {
        let a = vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)];
        let b = vec![Q::from_integer(4), Q::from_integer(5), Q::from_integer(6)];

        let result = parallel_dot(&a, &b);
        // 1*4 + 2*5 + 3*6 = 32
        assert_eq!(result, Q::from_integer(32));
    }

    #[test]
    fn test_parallel_add() {
        let a = vec![Q::from_integer(1), Q::from_integer(2)];
        let b = vec![Q::from_integer(3), Q::from_integer(4)];

        let c = parallel_add(&a, &b);
        assert_eq!(c, vec![Q::from_integer(4), Q::from_integer(6)]);
    }

    #[test]
    fn test_parallel_scale() {
        let alpha = Q::from_integer(2);
        let a = vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)];

        let b = parallel_scale(&alpha, &a);
        assert_eq!(
            b,
            vec![Q::from_integer(2), Q::from_integer(4), Q::from_integer(6)]
        );
    }
}
