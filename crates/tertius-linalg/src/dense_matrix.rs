//! Dense matrix implementation for small matrices.
//!
//! Dense matrices are more efficient for small sizes (< 100 rows)
//! due to better cache locality and simpler access patterns.

use std::ops::{Add, Index, IndexMut, Sub};

use num_traits::One;
use rayon::prelude::*;

use tertius_rings::traits::{Field, Ring};

/// Dense matrix stored in row-major order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseMatrix<R> {
    /// Matrix entries in row-major order.
    data: Vec<R>,
    /// Number of rows.
    num_rows: usize,
    /// Number of columns.
    num_cols: usize,
}

impl<R: Ring + Clone> DenseMatrix<R> {
    /// Creates a new matrix filled with zeros.
    #[must_use]
    pub fn zeros(num_rows: usize, num_cols: usize) -> Self {
        Self {
            data: vec![R::zero(); num_rows * num_cols],
            num_rows,
            num_cols,
        }
    }

    /// Creates a matrix from a 2D vector.
    #[must_use]
    pub fn from_rows(rows: Vec<Vec<R>>) -> Self {
        if rows.is_empty() {
            return Self::zeros(0, 0);
        }
        let num_rows = rows.len();
        let num_cols = rows[0].len();
        let data: Vec<R> = rows.into_iter().flatten().collect();
        assert_eq!(data.len(), num_rows * num_cols);
        Self {
            data,
            num_rows,
            num_cols,
        }
    }

    /// Creates an identity matrix.
    #[must_use]
    pub fn identity(n: usize) -> Self
    where
        R: One,
    {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = <R as One>::one();
        }
        m
    }

    /// Returns the number of rows.
    #[must_use]
    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    /// Returns the number of columns.
    #[must_use]
    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    /// Checks if the matrix is square.
    #[must_use]
    pub fn is_square(&self) -> bool {
        self.num_rows == self.num_cols
    }

    /// Returns a reference to the entry at (row, col).
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> Option<&R> {
        if row < self.num_rows && col < self.num_cols {
            Some(&self.data[row * self.num_cols + col])
        } else {
            None
        }
    }

    /// Returns a mutable reference to the entry at (row, col).
    pub fn get_mut(&mut self, row: usize, col: usize) -> Option<&mut R> {
        if row < self.num_rows && col < self.num_cols {
            Some(&mut self.data[row * self.num_cols + col])
        } else {
            None
        }
    }

    /// Returns a slice of the specified row.
    #[must_use]
    pub fn row(&self, row: usize) -> &[R] {
        let start = row * self.num_cols;
        &self.data[start..start + self.num_cols]
    }

    /// Returns a mutable slice of the specified row.
    pub fn row_mut(&mut self, row: usize) -> &mut [R] {
        let start = row * self.num_cols;
        &mut self.data[start..start + self.num_cols]
    }

    /// Returns a column as a vector.
    #[must_use]
    pub fn col(&self, col: usize) -> Vec<R> {
        (0..self.num_rows)
            .map(|row| self[(row, col)].clone())
            .collect()
    }

    /// Sets a column from a slice.
    pub fn set_col(&mut self, col: usize, values: &[R]) {
        assert_eq!(values.len(), self.num_rows);
        for (row, val) in values.iter().enumerate() {
            self[(row, col)] = val.clone();
        }
    }

    /// Matrix-vector multiply: y = A * x.
    #[must_use]
    pub fn mv(&self, x: &[R]) -> Vec<R> {
        assert_eq!(x.len(), self.num_cols);
        (0..self.num_rows)
            .map(|row| {
                self.row(row)
                    .iter()
                    .zip(x.iter())
                    .fold(R::zero(), |acc, (a, b)| acc + a.clone() * b.clone())
            })
            .collect()
    }

    /// Matrix-matrix multiply: C = A * B.
    #[must_use]
    pub fn mm(&self, other: &Self) -> Self {
        assert_eq!(self.num_cols, other.num_rows);

        let mut result = Self::zeros(self.num_rows, other.num_cols);
        for i in 0..self.num_rows {
            for j in 0..other.num_cols {
                let mut sum = R::zero();
                for k in 0..self.num_cols {
                    sum = sum + self[(i, k)].clone() * other[(k, j)].clone();
                }
                result[(i, j)] = sum;
            }
        }
        result
    }

    /// Returns the transpose of the matrix.
    #[must_use]
    pub fn transpose(&self) -> Self {
        let mut result = Self::zeros(self.num_cols, self.num_rows);
        for i in 0..self.num_rows {
            for j in 0..self.num_cols {
                result[(j, i)] = self[(i, j)].clone();
            }
        }
        result
    }

    /// Scales all entries by a scalar.
    #[must_use]
    pub fn scale(&self, scalar: &R) -> Self {
        Self {
            data: self.data.iter().map(|v| v.clone() * scalar.clone()).collect(),
            num_rows: self.num_rows,
            num_cols: self.num_cols,
        }
    }

    /// Swaps two rows in-place.
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }
        let i_start = i * self.num_cols;
        let j_start = j * self.num_cols;
        for k in 0..self.num_cols {
            self.data.swap(i_start + k, j_start + k);
        }
    }

    /// Adds a scaled row to another: row[target] += scale * row[source].
    pub fn add_scaled_row(&mut self, target: usize, source: usize, scale: &R) {
        for k in 0..self.num_cols {
            let val = self[(source, k)].clone() * scale.clone();
            self[(target, k)] = self[(target, k)].clone() + val;
        }
    }

    /// Scales a row by a scalar.
    pub fn scale_row(&mut self, row: usize, scale: &R) {
        for k in 0..self.num_cols {
            self[(row, k)] = self[(row, k)].clone() * scale.clone();
        }
    }
}

impl<R: Ring + Clone + Send + Sync> DenseMatrix<R> {
    /// Matrix-vector multiply (parallel): y = A * x.
    #[must_use]
    pub fn mv_parallel(&self, x: &[R]) -> Vec<R> {
        assert_eq!(x.len(), self.num_cols);
        (0..self.num_rows)
            .into_par_iter()
            .map(|row| {
                self.row(row)
                    .iter()
                    .zip(x.iter())
                    .fold(R::zero(), |acc, (a, b)| acc + a.clone() * b.clone())
            })
            .collect()
    }

    /// Matrix-matrix multiply (parallel): C = A * B.
    #[must_use]
    pub fn mm_parallel(&self, other: &Self) -> Self {
        assert_eq!(self.num_cols, other.num_rows);

        let data: Vec<R> = (0..self.num_rows)
            .into_par_iter()
            .flat_map(|i| {
                (0..other.num_cols)
                    .map(|j| {
                        let mut sum = R::zero();
                        for k in 0..self.num_cols {
                            sum = sum + self[(i, k)].clone() * other[(k, j)].clone();
                        }
                        sum
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        Self {
            data,
            num_rows: self.num_rows,
            num_cols: other.num_cols,
        }
    }
}

impl<R: Field + Clone + One> DenseMatrix<R> {
    /// Gaussian elimination with partial pivoting.
    ///
    /// Returns (row-echelon form, transformation matrix, rank).
    #[must_use]
    pub fn row_echelon(&self) -> (Self, Self, usize) {
        let mut m = self.clone();
        let mut transform = Self::identity(self.num_rows);
        let mut pivot_row = 0;
        let mut pivot_col = 0;

        while pivot_row < m.num_rows && pivot_col < m.num_cols {
            // Find pivot (first non-zero in column)
            let mut max_row = pivot_row;
            for row in pivot_row + 1..m.num_rows {
                if !m[(row, pivot_col)].is_zero() {
                    max_row = row;
                    break;
                }
            }

            if m[(max_row, pivot_col)].is_zero() {
                // No pivot in this column
                pivot_col += 1;
                continue;
            }

            // Swap rows if needed
            if max_row != pivot_row {
                m.swap_rows(pivot_row, max_row);
                transform.swap_rows(pivot_row, max_row);
            }

            // Scale pivot row to make pivot = 1
            let pivot_val = m[(pivot_row, pivot_col)].clone();
            if let Some(inv) = pivot_val.inv() {
                m.scale_row(pivot_row, &inv);
                transform.scale_row(pivot_row, &inv);
            }

            // Eliminate entries below pivot
            for row in pivot_row + 1..m.num_rows {
                if !m[(row, pivot_col)].is_zero() {
                    let factor = m[(row, pivot_col)].clone().neg();
                    m.add_scaled_row(row, pivot_row, &factor);
                    transform.add_scaled_row(row, pivot_row, &factor);
                }
            }

            pivot_row += 1;
            pivot_col += 1;
        }

        let rank = pivot_row;
        (m, transform, rank)
    }

    /// Reduced row echelon form (RREF) using Gauss-Jordan elimination.
    #[must_use]
    pub fn rref(&self) -> (Self, Self, usize) {
        let (mut m, mut transform, rank) = self.row_echelon();

        // Back-substitution to eliminate above pivots
        for pivot_row in (0..rank).rev() {
            // Find pivot column
            let mut pivot_col = 0;
            for col in 0..m.num_cols {
                if !m[(pivot_row, col)].is_zero() {
                    pivot_col = col;
                    break;
                }
            }

            // Eliminate entries above pivot
            for row in 0..pivot_row {
                if !m[(row, pivot_col)].is_zero() {
                    let factor = m[(row, pivot_col)].clone().neg();
                    m.add_scaled_row(row, pivot_row, &factor);
                    transform.add_scaled_row(row, pivot_row, &factor);
                }
            }
        }

        (m, transform, rank)
    }

    /// Computes the null space basis.
    ///
    /// Returns vectors that span the kernel of the matrix.
    #[must_use]
    pub fn null_space(&self) -> Vec<Vec<R>> {
        let (rref, _, rank) = self.rref();
        let mut basis = Vec::new();

        // Find pivot and free columns
        let mut pivot_cols = Vec::new();
        let mut row = 0;
        for col in 0..self.num_cols {
            if row < self.num_rows && !rref[(row, col)].is_zero() {
                pivot_cols.push(col);
                row += 1;
            }
        }

        // For each free variable, construct a null space vector
        for col in 0..self.num_cols {
            if pivot_cols.contains(&col) {
                continue;
            }

            let mut vec = vec![R::zero(); self.num_cols];
            vec[col] = <R as One>::one();

            // Fill in entries corresponding to pivot variables
            for (pivot_row, &pivot_col) in pivot_cols.iter().enumerate() {
                if pivot_col < self.num_cols {
                    vec[pivot_col] = rref[(pivot_row, col)].clone().neg();
                }
            }

            basis.push(vec);
        }

        basis
    }

    /// Solves the linear system Ax = b.
    ///
    /// Returns None if no solution exists.
    #[must_use]
    pub fn solve(&self, b: &[R]) -> Option<Vec<R>> {
        assert_eq!(b.len(), self.num_rows);

        // Augmented matrix [A | b]
        let mut aug = Self::zeros(self.num_rows, self.num_cols + 1);
        for i in 0..self.num_rows {
            for j in 0..self.num_cols {
                aug[(i, j)] = self[(i, j)].clone();
            }
            aug[(i, self.num_cols)] = b[i].clone();
        }

        let (rref, _, rank) = aug.rref();

        // Check for inconsistency: row of form [0 0 ... 0 | c] where c != 0
        for row in rank..self.num_rows {
            if !rref[(row, self.num_cols)].is_zero() {
                return None; // No solution
            }
        }

        // Extract solution (assuming unique solution for now)
        let mut x = vec![R::zero(); self.num_cols];
        for row in (0..rank).rev() {
            // Find pivot column
            let mut pivot_col = 0;
            for col in 0..self.num_cols {
                if !rref[(row, col)].is_zero() {
                    pivot_col = col;
                    break;
                }
            }
            x[pivot_col] = rref[(row, self.num_cols)].clone();
        }

        Some(x)
    }

    /// Computes the determinant (for square matrices).
    #[must_use]
    pub fn det(&self) -> R {
        assert!(self.is_square());
        let n = self.num_rows;

        if n == 0 {
            return <R as One>::one();
        }
        if n == 1 {
            return self[(0, 0)].clone();
        }
        if n == 2 {
            return self[(0, 0)].clone() * self[(1, 1)].clone()
                + (self[(0, 1)].clone() * self[(1, 0)].clone()).neg();
        }

        // LU decomposition approach
        let mut m = self.clone();
        let mut det = <R as One>::one();
        let mut sign = <R as One>::one();

        for col in 0..n {
            // Find pivot
            let mut pivot_row = col;
            while pivot_row < n && m[(pivot_row, col)].is_zero() {
                pivot_row += 1;
            }

            if pivot_row == n {
                return R::zero(); // Singular
            }

            if pivot_row != col {
                m.swap_rows(col, pivot_row);
                sign = sign.neg();
            }

            let pivot = m[(col, col)].clone();
            det = det * pivot.clone();

            // Eliminate below
            if let Some(inv) = pivot.inv() {
                for row in col + 1..n {
                    if !m[(row, col)].is_zero() {
                        let factor = m[(row, col)].clone() * inv.clone();
                        for k in col..n {
                            let val = m[(col, k)].clone() * factor.clone();
                            m[(row, k)] = m[(row, k)].clone() + val.neg();
                        }
                    }
                }
            }
        }

        det * sign
    }

    /// Computes the inverse (for square matrices).
    #[must_use]
    pub fn inverse(&self) -> Option<Self> {
        assert!(self.is_square());
        let n = self.num_rows;

        // Augmented matrix [A | I]
        let mut aug = Self::zeros(n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                aug[(i, j)] = self[(i, j)].clone();
            }
            aug[(i, n + i)] = <R as One>::one();
        }

        let (rref, _, rank) = aug.rref();

        if rank != n {
            return None; // Singular
        }

        // Extract inverse from right half
        let mut inv = Self::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                inv[(i, j)] = rref[(i, n + j)].clone();
            }
        }

        Some(inv)
    }
}

impl<R> Index<(usize, usize)> for DenseMatrix<R> {
    type Output = R;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[row * self.num_cols + col]
    }
}

impl<R> IndexMut<(usize, usize)> for DenseMatrix<R> {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[row * self.num_cols + col]
    }
}

impl<R: Ring + Clone + Add<Output = R>> Add for &DenseMatrix<R> {
    type Output = DenseMatrix<R>;

    fn add(self, other: Self) -> DenseMatrix<R> {
        assert_eq!(self.num_rows, other.num_rows);
        assert_eq!(self.num_cols, other.num_cols);

        DenseMatrix {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(a, b)| a.clone() + b.clone())
                .collect(),
            num_rows: self.num_rows,
            num_cols: self.num_cols,
        }
    }
}

impl<R: Ring + Clone + Sub<Output = R>> Sub for &DenseMatrix<R> {
    type Output = DenseMatrix<R>;

    fn sub(self, other: Self) -> DenseMatrix<R> {
        assert_eq!(self.num_rows, other.num_rows);
        assert_eq!(self.num_cols, other.num_cols);

        DenseMatrix {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(a, b)| a.clone() - b.clone())
                .collect(),
            num_rows: self.num_rows,
            num_cols: self.num_cols,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::integers::Z;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_zeros() {
        let m: DenseMatrix<Z> = DenseMatrix::zeros(3, 4);
        assert_eq!(m.num_rows(), 3);
        assert_eq!(m.num_cols(), 4);
        for i in 0..3 {
            for j in 0..4 {
                assert_eq!(m[(i, j)], Z::new(0));
            }
        }
    }

    #[test]
    fn test_identity() {
        let id: DenseMatrix<Z> = DenseMatrix::identity(3);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(id[(i, j)], Z::new(1));
                } else {
                    assert_eq!(id[(i, j)], Z::new(0));
                }
            }
        }
    }

    #[test]
    fn test_mv() {
        let m = DenseMatrix::from_rows(vec![
            vec![Z::new(1), Z::new(2), Z::new(3)],
            vec![Z::new(4), Z::new(5), Z::new(6)],
        ]);
        let x = vec![Z::new(1), Z::new(2), Z::new(3)];
        let y = m.mv(&x);
        // [1*1 + 2*2 + 3*3, 4*1 + 5*2 + 6*3] = [14, 32]
        assert_eq!(y, vec![Z::new(14), Z::new(32)]);
    }

    #[test]
    fn test_mm() {
        let a = DenseMatrix::from_rows(vec![
            vec![Z::new(1), Z::new(2)],
            vec![Z::new(3), Z::new(4)],
        ]);
        let b = DenseMatrix::from_rows(vec![
            vec![Z::new(5), Z::new(6)],
            vec![Z::new(7), Z::new(8)],
        ]);
        let c = a.mm(&b);
        // [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]] = [[19, 22], [43, 50]]
        assert_eq!(c[(0, 0)], Z::new(19));
        assert_eq!(c[(0, 1)], Z::new(22));
        assert_eq!(c[(1, 0)], Z::new(43));
        assert_eq!(c[(1, 1)], Z::new(50));
    }

    #[test]
    fn test_transpose() {
        let m = DenseMatrix::from_rows(vec![
            vec![Z::new(1), Z::new(2), Z::new(3)],
            vec![Z::new(4), Z::new(5), Z::new(6)],
        ]);
        let t = m.transpose();
        assert_eq!(t.num_rows(), 3);
        assert_eq!(t.num_cols(), 2);
        assert_eq!(t[(0, 0)], Z::new(1));
        assert_eq!(t[(1, 0)], Z::new(2));
        assert_eq!(t[(2, 1)], Z::new(6));
    }

    #[test]
    fn test_det() {
        // 2x2 determinant
        let m = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(3), Q::from_integer(8)],
            vec![Q::from_integer(4), Q::from_integer(6)],
        ]);
        let det = m.det();
        // 3*6 - 8*4 = 18 - 32 = -14
        assert_eq!(det, Q::from_integer(-14));
    }

    #[test]
    fn test_inverse() {
        let m = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(4), Q::from_integer(7)],
            vec![Q::from_integer(2), Q::from_integer(6)],
        ]);
        let inv = m.inverse().unwrap();
        let product = m.mm(&inv);

        // Should be identity
        assert_eq!(product[(0, 0)], Q::from_integer(1));
        assert_eq!(product[(0, 1)], Q::from_integer(0));
        assert_eq!(product[(1, 0)], Q::from_integer(0));
        assert_eq!(product[(1, 1)], Q::from_integer(1));
    }

    #[test]
    fn test_solve() {
        // Solve Ax = b where A = [[1, 2], [3, 4]], b = [5, 11]
        // Solution: x = [1, 2]
        let a = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(2)],
            vec![Q::from_integer(3), Q::from_integer(4)],
        ]);
        let b = vec![Q::from_integer(5), Q::from_integer(11)];
        let x = a.solve(&b).unwrap();

        assert_eq!(x[0], Q::from_integer(1));
        assert_eq!(x[1], Q::from_integer(2));
    }

    #[test]
    fn test_null_space() {
        // Matrix with non-trivial null space
        let m = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)],
            vec![Q::from_integer(2), Q::from_integer(4), Q::from_integer(6)],
        ]);
        let null = m.null_space();

        // Should have 2 basis vectors (rank 1, 3 cols - 1 = 2 free vars)
        // Actually rank is 1, so null space dimension is 3 - 1 = 2
        assert!(!null.is_empty());

        // Verify each vector is in null space
        for v in &null {
            let result = m.mv(v);
            for r in result {
                assert!(r.is_zero());
            }
        }
    }
}
