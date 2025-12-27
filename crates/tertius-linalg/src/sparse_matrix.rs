//! Sparse matrix in Compressed Sparse Row (CSR) format.
//!
//! CSR format is optimal for row-wise access patterns, which are common in:
//! - Sparse matrix-vector multiplication (SpMV)
//! - Gaussian elimination (row operations)
//! - Iterative solvers (Block Wiedemann)

use std::ops::{Add, Neg, Sub};

use num_traits::One;
use rayon::prelude::*;

use tertius_rings::traits::Ring;

/// Sparse matrix in Compressed Sparse Row (CSR) format.
///
/// # Memory Layout
///
/// For an m×n matrix with nnz non-zero entries:
/// - `values`: Vec of nnz non-zero values
/// - `col_indices`: Vec of nnz column indices
/// - `row_ptrs`: Vec of m+1 row pointers
///
/// Row i contains entries from `row_ptrs[i]` to `row_ptrs[i+1]`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CsrMatrix<R> {
    /// Non-zero values in row-major order.
    values: Vec<R>,
    /// Column index for each non-zero value.
    col_indices: Vec<usize>,
    /// Row pointers: row_ptrs[i] is the index into values where row i starts.
    row_ptrs: Vec<usize>,
    /// Number of columns.
    num_cols: usize,
}

impl<R: Ring + Clone> CsrMatrix<R> {
    /// Creates a new empty sparse matrix.
    #[must_use]
    pub fn new(num_rows: usize, num_cols: usize) -> Self {
        Self {
            values: Vec::new(),
            col_indices: Vec::new(),
            row_ptrs: vec![0; num_rows + 1],
            num_cols,
        }
    }

    /// Creates a sparse matrix from a dense matrix.
    ///
    /// Zero entries are not stored.
    #[must_use]
    pub fn from_dense(dense: &[Vec<R>]) -> Self {
        if dense.is_empty() {
            return Self::new(0, 0);
        }

        let num_rows = dense.len();
        let num_cols = dense[0].len();
        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptrs = Vec::with_capacity(num_rows + 1);

        for row in dense {
            row_ptrs.push(values.len());
            for (col, val) in row.iter().enumerate() {
                if !val.is_zero() {
                    values.push(val.clone());
                    col_indices.push(col);
                }
            }
        }
        row_ptrs.push(values.len());

        Self {
            values,
            col_indices,
            row_ptrs,
            num_cols,
        }
    }

    /// Creates a sparse matrix from triplets (row, col, value).
    ///
    /// Duplicate entries are summed. Entries are sorted by (row, col).
    #[must_use]
    pub fn from_triplets(
        num_rows: usize,
        num_cols: usize,
        triplets: &[(usize, usize, R)],
    ) -> Self {
        if triplets.is_empty() {
            return Self::new(num_rows, num_cols);
        }

        // Sort triplets by (row, col)
        let mut sorted: Vec<_> = triplets.to_vec();
        sorted.sort_by_key(|(r, c, _)| (*r, *c));

        // Build CSR directly from sorted triplets
        let mut values: Vec<R> = Vec::new();
        let mut col_indices: Vec<usize> = Vec::new();
        let mut row_ptrs = Vec::with_capacity(num_rows + 1);
        row_ptrs.push(0);

        let mut prev_row = usize::MAX;
        let mut prev_col = usize::MAX;
        let mut current_row = 0;

        for (row, col, val) in sorted {
            if row == prev_row && col == prev_col {
                // Sum duplicate entries
                if let Some(last) = values.last_mut() {
                    *last = last.clone() + val;
                }
            } else {
                // Fill in empty rows
                while current_row < row {
                    row_ptrs.push(values.len());
                    current_row += 1;
                }

                if !val.is_zero() {
                    values.push(val);
                    col_indices.push(col);
                }
                prev_row = row;
                prev_col = col;
            }
        }

        // Fill remaining row pointers
        while row_ptrs.len() <= num_rows {
            row_ptrs.push(values.len());
        }

        Self {
            values,
            col_indices,
            row_ptrs,
            num_cols,
        }
    }

    /// Creates an identity matrix of size n×n.
    #[must_use]
    pub fn identity(n: usize) -> Self
    where
        R: One,
    {
        let values: Vec<R> = (0..n).map(|_| <R as One>::one()).collect();
        let col_indices: Vec<usize> = (0..n).collect();
        let row_ptrs: Vec<usize> = (0..=n).collect();

        Self {
            values,
            col_indices,
            row_ptrs,
            num_cols: n,
        }
    }

    /// Returns the number of rows.
    #[must_use]
    pub fn num_rows(&self) -> usize {
        self.row_ptrs.len().saturating_sub(1)
    }

    /// Returns the number of columns.
    #[must_use]
    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    /// Returns the number of non-zero entries.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Returns the density (fraction of non-zero entries).
    #[must_use]
    pub fn density(&self) -> f64 {
        let total = self.num_rows() * self.num_cols;
        if total == 0 {
            0.0
        } else {
            self.nnz() as f64 / total as f64
        }
    }

    /// Checks if the matrix is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.num_rows() == 0 || self.num_cols == 0
    }

    /// Returns an iterator over non-zero entries in a row.
    pub fn row_iter(&self, row: usize) -> impl Iterator<Item = (usize, &R)> {
        let start = self.row_ptrs[row];
        let end = self.row_ptrs[row + 1];
        self.col_indices[start..end]
            .iter()
            .zip(&self.values[start..end])
            .map(|(&col, val)| (col, val))
    }

    /// Returns the entry at (row, col), or None if zero.
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> Option<&R> {
        let start = self.row_ptrs[row];
        let end = self.row_ptrs[row + 1];

        // Binary search for the column
        let idx = self.col_indices[start..end]
            .binary_search(&col)
            .ok()?;
        Some(&self.values[start + idx])
    }

    /// Computes the dot product of row i with a vector.
    fn row_dot(&self, row: usize, x: &[R]) -> R {
        let start = self.row_ptrs[row];
        let end = self.row_ptrs[row + 1];

        let mut result = R::zero();
        for (&col, val) in self.col_indices[start..end]
            .iter()
            .zip(&self.values[start..end])
        {
            result = result + val.clone() * x[col].clone();
        }
        result
    }

    /// Sparse matrix-vector multiply: y = A * x (sequential).
    #[must_use]
    pub fn spmv(&self, x: &[R]) -> Vec<R> {
        assert_eq!(x.len(), self.num_cols, "Vector dimension mismatch");

        (0..self.num_rows())
            .map(|i| self.row_dot(i, x))
            .collect()
    }
}

impl<R: Ring + Clone + Send + Sync> CsrMatrix<R> {
    /// Sparse matrix-vector multiply: y = A * x (parallel).
    ///
    /// Uses rayon for parallel row computation.
    #[must_use]
    pub fn spmv_parallel(&self, x: &[R]) -> Vec<R> {
        assert_eq!(x.len(), self.num_cols, "Vector dimension mismatch");

        (0..self.num_rows())
            .into_par_iter()
            .map(|i| self.row_dot(i, x))
            .collect()
    }

    /// Computes A^T * x (transpose-vector multiply) in parallel.
    #[must_use]
    pub fn spmv_transpose_parallel(&self, x: &[R]) -> Vec<R> {
        assert_eq!(x.len(), self.num_rows(), "Vector dimension mismatch");

        // Use atomic-like accumulation via parallel fold
        let partial_results: Vec<Vec<(usize, R)>> = (0..self.num_rows())
            .into_par_iter()
            .filter_map(|i| {
                if x[i].is_zero() {
                    return None;
                }
                let contributions: Vec<_> = self
                    .row_iter(i)
                    .map(|(col, val)| (col, val.clone() * x[i].clone()))
                    .collect();
                if contributions.is_empty() {
                    None
                } else {
                    Some(contributions)
                }
            })
            .collect();

        // Aggregate results
        let mut result = vec![R::zero(); self.num_cols];
        for contributions in partial_results {
            for (col, val) in contributions {
                result[col] = result[col].clone() + val;
            }
        }
        result
    }
}

impl<R: Ring + Clone> CsrMatrix<R> {
    /// Converts to dense matrix representation.
    #[must_use]
    pub fn to_dense(&self) -> Vec<Vec<R>> {
        let mut dense = vec![vec![R::zero(); self.num_cols]; self.num_rows()];
        for row in 0..self.num_rows() {
            for (col, val) in self.row_iter(row) {
                dense[row][col] = val.clone();
            }
        }
        dense
    }

    /// Returns the transpose of the matrix.
    #[must_use]
    pub fn transpose(&self) -> Self {
        // Build CSR of transpose via triplets
        let triplets: Vec<_> = (0..self.num_rows())
            .flat_map(|row| {
                self.row_iter(row)
                    .map(move |(col, val)| (col, row, val.clone()))
            })
            .collect();

        Self::from_triplets(self.num_cols, self.num_rows(), &triplets)
    }

    /// Scales all entries by a scalar.
    #[must_use]
    pub fn scale(&self, scalar: &R) -> Self {
        Self {
            values: self.values.iter().map(|v| v.clone() * scalar.clone()).collect(),
            col_indices: self.col_indices.clone(),
            row_ptrs: self.row_ptrs.clone(),
            num_cols: self.num_cols,
        }
    }
}

impl<R: Ring + Clone + Add<Output = R>> Add for &CsrMatrix<R> {
    type Output = CsrMatrix<R>;

    fn add(self, other: Self) -> CsrMatrix<R> {
        assert_eq!(self.num_rows(), other.num_rows());
        assert_eq!(self.num_cols, other.num_cols);

        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptrs = Vec::with_capacity(self.num_rows() + 1);

        for row in 0..self.num_rows() {
            row_ptrs.push(values.len());

            let mut i = self.row_ptrs[row];
            let mut j = other.row_ptrs[row];
            let i_end = self.row_ptrs[row + 1];
            let j_end = other.row_ptrs[row + 1];

            while i < i_end && j < j_end {
                let col_i = self.col_indices[i];
                let col_j = other.col_indices[j];

                match col_i.cmp(&col_j) {
                    std::cmp::Ordering::Less => {
                        values.push(self.values[i].clone());
                        col_indices.push(col_i);
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        values.push(other.values[j].clone());
                        col_indices.push(col_j);
                        j += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        let sum = self.values[i].clone() + other.values[j].clone();
                        if !sum.is_zero() {
                            values.push(sum);
                            col_indices.push(col_i);
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }

            // Append remaining entries
            while i < i_end {
                values.push(self.values[i].clone());
                col_indices.push(self.col_indices[i]);
                i += 1;
            }
            while j < j_end {
                values.push(other.values[j].clone());
                col_indices.push(other.col_indices[j]);
                j += 1;
            }
        }
        row_ptrs.push(values.len());

        CsrMatrix {
            values,
            col_indices,
            row_ptrs,
            num_cols: self.num_cols,
        }
    }
}

impl<R: Ring + Clone + Sub<Output = R> + Neg<Output = R>> Sub for &CsrMatrix<R> {
    type Output = CsrMatrix<R>;

    fn sub(self, other: Self) -> CsrMatrix<R> {
        // A - B = A + (-B)
        let neg_other = other.scale(&R::zero().neg().neg()); // This doesn't work, need actual neg
        self + &CsrMatrix {
            values: other.values.iter().map(|v| v.clone().neg()).collect(),
            col_indices: other.col_indices.clone(),
            row_ptrs: other.row_ptrs.clone(),
            num_cols: other.num_cols,
        }
    }
}

/// Iterator over rows of a CSR matrix.
pub struct RowIter<'a, R> {
    matrix: &'a CsrMatrix<R>,
    current_row: usize,
}

impl<'a, R: Ring + Clone> Iterator for RowIter<'a, R> {
    type Item = Vec<(usize, R)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_row >= self.matrix.num_rows() {
            return None;
        }

        let row: Vec<_> = self
            .matrix
            .row_iter(self.current_row)
            .map(|(col, val)| (col, val.clone()))
            .collect();
        self.current_row += 1;
        Some(row)
    }
}

impl<R: Ring + Clone> CsrMatrix<R> {
    /// Returns an iterator over all rows.
    #[must_use]
    pub fn rows(&self) -> RowIter<'_, R> {
        RowIter {
            matrix: self,
            current_row: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::integers::Z;

    #[test]
    fn test_from_dense() {
        let dense = vec![
            vec![Z::new(1), Z::new(0), Z::new(2)],
            vec![Z::new(0), Z::new(3), Z::new(0)],
            vec![Z::new(4), Z::new(0), Z::new(5)],
        ];
        let sparse = CsrMatrix::from_dense(&dense);

        assert_eq!(sparse.num_rows(), 3);
        assert_eq!(sparse.num_cols(), 3);
        assert_eq!(sparse.nnz(), 5);

        assert_eq!(sparse.get(0, 0), Some(&Z::new(1)));
        assert_eq!(sparse.get(0, 1), None);
        assert_eq!(sparse.get(0, 2), Some(&Z::new(2)));
        assert_eq!(sparse.get(1, 1), Some(&Z::new(3)));
        assert_eq!(sparse.get(2, 0), Some(&Z::new(4)));
        assert_eq!(sparse.get(2, 2), Some(&Z::new(5)));
    }

    #[test]
    fn test_spmv() {
        let dense = vec![
            vec![Z::new(1), Z::new(2), Z::new(0)],
            vec![Z::new(0), Z::new(3), Z::new(4)],
            vec![Z::new(5), Z::new(0), Z::new(6)],
        ];
        let sparse = CsrMatrix::from_dense(&dense);
        let x = vec![Z::new(1), Z::new(2), Z::new(3)];

        let result = sparse.spmv(&x);
        // [1*1 + 2*2 + 0*3, 0*1 + 3*2 + 4*3, 5*1 + 0*2 + 6*3]
        // = [5, 18, 23]
        assert_eq!(result, vec![Z::new(5), Z::new(18), Z::new(23)]);
    }

    #[test]
    fn test_spmv_parallel() {
        let dense = vec![
            vec![Z::new(1), Z::new(2), Z::new(0)],
            vec![Z::new(0), Z::new(3), Z::new(4)],
            vec![Z::new(5), Z::new(0), Z::new(6)],
        ];
        let sparse = CsrMatrix::from_dense(&dense);
        let x = vec![Z::new(1), Z::new(2), Z::new(3)];

        let result = sparse.spmv_parallel(&x);
        assert_eq!(result, vec![Z::new(5), Z::new(18), Z::new(23)]);
    }

    #[test]
    fn test_transpose() {
        let dense = vec![
            vec![Z::new(1), Z::new(2), Z::new(0)],
            vec![Z::new(0), Z::new(3), Z::new(4)],
        ];
        let sparse = CsrMatrix::from_dense(&dense);
        let transposed = sparse.transpose();

        assert_eq!(transposed.num_rows(), 3);
        assert_eq!(transposed.num_cols(), 2);
        assert_eq!(transposed.get(0, 0), Some(&Z::new(1)));
        assert_eq!(transposed.get(1, 0), Some(&Z::new(2)));
        assert_eq!(transposed.get(1, 1), Some(&Z::new(3)));
        assert_eq!(transposed.get(2, 1), Some(&Z::new(4)));
    }

    #[test]
    fn test_identity() {
        let id: CsrMatrix<Z> = CsrMatrix::identity(3);
        let x = vec![Z::new(1), Z::new(2), Z::new(3)];
        let result = id.spmv(&x);
        assert_eq!(result, x);
    }

    #[test]
    fn test_add() {
        let a = CsrMatrix::from_dense(&vec![
            vec![Z::new(1), Z::new(0)],
            vec![Z::new(0), Z::new(2)],
        ]);
        let b = CsrMatrix::from_dense(&vec![
            vec![Z::new(0), Z::new(3)],
            vec![Z::new(4), Z::new(0)],
        ]);
        let c = &a + &b;

        assert_eq!(c.get(0, 0), Some(&Z::new(1)));
        assert_eq!(c.get(0, 1), Some(&Z::new(3)));
        assert_eq!(c.get(1, 0), Some(&Z::new(4)));
        assert_eq!(c.get(1, 1), Some(&Z::new(2)));
    }

    #[test]
    fn test_to_dense() {
        let original = vec![
            vec![Z::new(1), Z::new(0), Z::new(2)],
            vec![Z::new(0), Z::new(3), Z::new(0)],
            vec![Z::new(4), Z::new(0), Z::new(5)],
        ];
        let sparse = CsrMatrix::from_dense(&original);
        let reconstructed = sparse.to_dense();
        assert_eq!(original, reconstructed);
    }
}
