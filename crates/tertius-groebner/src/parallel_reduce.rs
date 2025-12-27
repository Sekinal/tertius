//! Parallel row reduction for Macaulay matrices.
//!
//! This module implements tiled parallel Gaussian elimination
//! optimized for sparse polynomial matrices.

use crate::macaulay::{MacaulayMatrix, MacaulayRow};
use num_traits::Zero;
use rayon::prelude::*;
use std::ops::Neg;
use tertius_rings::traits::{Field, Ring};

/// Result of reducing a Macaulay matrix.
pub struct ReductionResult<R> {
    /// Reduced rows (in echelon form).
    pub reduced_rows: Vec<ReducedRow<R>>,
    /// Number of reduction steps performed.
    pub reduction_count: usize,
}

/// A reduced row with its pivot column.
#[derive(Clone, Debug)]
pub struct ReducedRow<R> {
    /// Pivot column (leading column).
    pub pivot: usize,
    /// Non-zero entries after pivot: (column, coefficient).
    /// Does not include the pivot itself (assumed to be 1 for monic).
    pub tail: Vec<(usize, R)>,
    /// Original source index.
    pub source_index: usize,
}

impl<R: Clone> ReducedRow<R> {
    /// Creates a new reduced row.
    pub fn new(pivot: usize, tail: Vec<(usize, R)>, source_index: usize) -> Self {
        Self {
            pivot,
            tail,
            source_index,
        }
    }

    /// Returns all entries including pivot with coefficient 1.
    pub fn to_entries(&self) -> Vec<(usize, R)>
    where
        R: num_traits::One + Clone,
    {
        let mut entries = vec![(self.pivot, R::one())];
        entries.extend(self.tail.iter().cloned());
        entries
    }
}

/// Reduces a Macaulay matrix to row echelon form using parallel row operations.
///
/// This is a sparse implementation that only stores non-zero entries.
pub fn reduce_matrix<R>(matrix: &MacaulayMatrix<R>) -> ReductionResult<R>
where
    R: Field + Clone + Send + Sync + Neg<Output = R>,
{
    let num_rows = matrix.num_rows();
    let num_cols = matrix.num_cols();

    if num_rows == 0 || num_cols == 0 {
        return ReductionResult {
            reduced_rows: vec![],
            reduction_count: 0,
        };
    }

    // Convert rows to mutable sparse format
    let mut rows: Vec<Option<SparseRow<R>>> = matrix
        .rows
        .iter()
        .map(|r| Some(SparseRow::from_macaulay(r)))
        .collect();

    let mut reduced_rows: Vec<ReducedRow<R>> = Vec::new();
    let mut reduction_count = 0;

    // Process each column
    for pivot_col in 0..num_cols {
        // Find a row with this pivot
        let mut pivot_row_idx = None;
        for (i, row) in rows.iter().enumerate() {
            if let Some(ref r) = row {
                if r.leading_column() == Some(pivot_col) {
                    pivot_row_idx = Some(i);
                    break;
                }
            }
        }

        let pivot_row_idx = match pivot_row_idx {
            Some(i) => i,
            None => continue,
        };

        // Extract and normalize the pivot row
        let pivot_row = rows[pivot_row_idx].take().unwrap();
        let pivot_row = pivot_row.make_monic();

        // Reduce other rows with this pivot (in parallel)
        let pivot_col_val = pivot_col;
        let pivot_tail: Vec<_> = pivot_row.entries[1..].to_vec();

        let reductions: usize = rows
            .par_iter_mut()
            .map(|row_opt| {
                if let Some(ref mut row) = row_opt {
                    if let Some(coeff) = row.coefficient_at(pivot_col_val) {
                        // Subtract coeff * pivot_row from this row
                        row.subtract_multiple(pivot_col_val, &pivot_tail, &coeff);
                        return 1;
                    }
                }
                0
            })
            .sum();
        reduction_count += reductions;

        // Store the reduced row
        reduced_rows.push(ReducedRow::new(
            pivot_col,
            pivot_tail.to_vec(),
            pivot_row.source_index,
        ));
    }

    ReductionResult {
        reduced_rows,
        reduction_count,
    }
}

/// A mutable sparse row for reduction.
#[derive(Clone)]
struct SparseRow<R> {
    /// Non-zero entries sorted by column.
    entries: Vec<(usize, R)>,
    /// Source polynomial index.
    source_index: usize,
}

impl<R: Field + Clone + Neg<Output = R>> SparseRow<R> {
    fn from_macaulay(row: &MacaulayRow<R>) -> Self {
        Self {
            entries: row.entries().to_vec(),
            source_index: row.source_index,
        }
    }

    fn leading_column(&self) -> Option<usize> {
        self.entries.first().map(|(col, _)| *col)
    }

    fn coefficient_at(&self, col: usize) -> Option<R> {
        self.entries
            .iter()
            .find(|(c, _)| *c == col)
            .map(|(_, coeff)| coeff.clone())
    }

    fn make_monic(mut self) -> Self {
        if let Some((_, lead_coeff)) = self.entries.first() {
            let inv = lead_coeff.inv().expect("leading coefficient must be invertible");
            for (_, coeff) in &mut self.entries {
                *coeff = coeff.clone() * inv.clone();
            }
        }
        self
    }

    fn subtract_multiple(&mut self, pivot_col: usize, pivot_tail: &[(usize, R)], multiplier: &R) {
        // Remove the entry at pivot_col
        self.entries.retain(|(c, _)| *c != pivot_col);

        // Subtract multiplier * each pivot_tail entry
        let neg_mult = multiplier.clone().neg();
        for (col, coeff) in pivot_tail {
            let adjustment = coeff.clone() * neg_mult.clone();
            self.add_to_column(*col, adjustment);
        }

        // Remove zeros
        self.entries.retain(|(_, c)| !c.is_zero());

        // Re-sort by column
        self.entries.sort_by_key(|(c, _)| *c);
    }

    fn add_to_column(&mut self, col: usize, value: R) {
        // Find existing entry
        for (c, coeff) in &mut self.entries {
            if *c == col {
                *coeff = coeff.clone() + value;
                return;
            }
        }
        // Add new entry
        self.entries.push((col, value));
    }
}

/// Interreduces a set of polynomials (makes each tail reduced w.r.t. others).
pub fn interreduce<R>(rows: &mut [ReducedRow<R>])
where
    R: Field + Clone + Neg<Output = R>,
{
    // Sort by pivot column
    rows.sort_by_key(|r| r.pivot);

    // For each row, reduce its tail by all previous pivots
    for i in 0..rows.len() {
        for j in 0..i {
            let pivot_col = rows[j].pivot;
            // Clone row j's tail to avoid borrow conflicts
            let row_j_tail: Vec<_> = rows[j].tail.clone();

            // Check if row i has this column in its tail
            if let Some(pos) = rows[i].tail.iter().position(|(c, _)| *c == pivot_col) {
                let coeff = rows[i].tail[pos].1.clone();

                // Remove the entry
                rows[i].tail.remove(pos);

                // Subtract coeff * row_j.tail from row_i.tail
                let neg_coeff = coeff.neg();
                for (col, c) in &row_j_tail {
                    let adjustment = c.clone() * neg_coeff.clone();

                    if let Some(entry) = rows[i].tail.iter_mut().find(|(cc, _)| cc == col) {
                        entry.1 = entry.1.clone() + adjustment;
                    } else {
                        rows[i].tail.push((*col, adjustment));
                    }
                }

                // Remove zeros and re-sort
                rows[i].tail.retain(|(_, c)| !c.is_zero());
                rows[i].tail.sort_by_key(|(c, _)| *c);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::finite_field::FiniteField;

    type GF7 = FiniteField<7>;

    #[test]
    fn test_sparse_row_operations() {
        let row = SparseRow {
            entries: vec![(0, GF7::new(2)), (2, GF7::new(3))],
            source_index: 0,
        };

        let monic = row.make_monic();
        assert_eq!(monic.entries[0].1, GF7::new(1));
        // 3 * inv(2) mod 7 = 3 * 4 mod 7 = 12 mod 7 = 5
        assert_eq!(monic.entries[1].1, GF7::new(5));
    }
}
