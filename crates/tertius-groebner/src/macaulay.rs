//! Macaulay matrix construction for F4-style reduction.
//!
//! This module builds sparse matrices from polynomials and their multiples,
//! enabling batch reduction via linear algebra.

use crate::labeled_poly::LabeledPoly;
use crate::monomial::PackedMonomial;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use tertius_linalg::sparse_matrix::CsrMatrix;
use tertius_rings::traits::Ring;

/// A row in the Macaulay matrix representing a polynomial.
#[derive(Clone, Debug)]
pub struct MacaulayRow<R> {
    /// Non-zero entries: (column index, coefficient).
    entries: Vec<(usize, R)>,
    /// Source polynomial index.
    pub source_index: usize,
    /// Multiplier monomial applied to the source.
    pub multiplier: PackedMonomial,
}

impl<R: Clone> MacaulayRow<R> {
    /// Creates a new row.
    pub fn new(entries: Vec<(usize, R)>, source_index: usize, multiplier: PackedMonomial) -> Self {
        Self {
            entries,
            source_index,
            multiplier,
        }
    }

    /// Returns the entries.
    pub fn entries(&self) -> &[(usize, R)] {
        &self.entries
    }

    /// Returns the leading column (smallest index with non-zero entry).
    pub fn leading_column(&self) -> Option<usize> {
        self.entries.first().map(|(col, _)| *col)
    }
}

/// The Macaulay matrix for a batch of S-polynomials.
pub struct MacaulayMatrix<R> {
    /// Rows of the matrix (polynomial multiples).
    pub rows: Vec<MacaulayRow<R>>,
    /// Column ordering: column index -> monomial.
    pub columns: Vec<PackedMonomial>,
    /// Inverse map: monomial -> column index.
    pub monomial_to_col: FxHashMap<MonomialKey, usize>,
}

/// Hashable key for monomials.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct MonomialKey(pub [u16; 16], pub u32);

impl From<&PackedMonomial> for MonomialKey {
    fn from(m: &PackedMonomial) -> Self {
        let mut arr = [0u16; 16];
        let exps = m.exponents();
        let n = exps.len().min(16);
        arr[..n].copy_from_slice(&exps[..n]);
        Self(arr, m.total_degree())
    }
}

impl<R: Ring + Clone + Send + Sync> MacaulayMatrix<R> {
    /// Constructs a Macaulay matrix from a set of labeled polynomials and required multiples.
    ///
    /// # Arguments
    /// - `polys`: The polynomials to include.
    /// - `multiples`: For each polynomial, the monomials to multiply by.
    pub fn construct(
        polys: &[LabeledPoly<R>],
        multiples: &[Vec<PackedMonomial>],
    ) -> Self {
        // Collect all monomials that will appear
        let mut all_monomials: Vec<PackedMonomial> = Vec::new();
        let mut seen: FxHashMap<MonomialKey, ()> = FxHashMap::default();

        for (poly, mults) in polys.iter().zip(multiples.iter()) {
            for mult in mults {
                for (_, m) in poly.terms() {
                    let product = m.mul(mult);
                    let key = MonomialKey::from(&product);
                    if !seen.contains_key(&key) {
                        seen.insert(key, ());
                        all_monomials.push(product);
                    }
                }
            }
        }

        // Sort monomials by grevlex descending (highest first = column 0)
        all_monomials.sort_by(|a, b| b.cmp_grevlex(a));

        // Build column index map
        let mut monomial_to_col: FxHashMap<MonomialKey, usize> = FxHashMap::default();
        for (col, m) in all_monomials.iter().enumerate() {
            monomial_to_col.insert(MonomialKey::from(m), col);
        }

        // Build rows in parallel
        let rows: Vec<MacaulayRow<R>> = polys
            .par_iter()
            .enumerate()
            .flat_map(|(poly_idx, poly)| {
                multiples[poly_idx]
                    .par_iter()
                    .map(|mult| {
                        let mut entries: Vec<(usize, R)> = Vec::new();

                        for (coeff, m) in poly.terms() {
                            let product = m.mul(mult);
                            let key = MonomialKey::from(&product);
                            if let Some(&col) = monomial_to_col.get(&key) {
                                entries.push((col, coeff.clone()));
                            }
                        }

                        // Sort by column
                        entries.sort_by_key(|(col, _)| *col);

                        MacaulayRow::new(entries, poly_idx, mult.clone())
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        Self {
            rows,
            columns: all_monomials,
            monomial_to_col,
        }
    }

    /// Converts to a CSR sparse matrix.
    pub fn to_csr(&self) -> CsrMatrix<R> {
        let num_rows = self.rows.len();
        let num_cols = self.columns.len();

        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptrs = Vec::with_capacity(num_rows + 1);
        row_ptrs.push(0);

        for row in &self.rows {
            for (col, val) in &row.entries {
                col_indices.push(*col);
                values.push(val.clone());
            }
            row_ptrs.push(values.len());
        }

        CsrMatrix::from_raw(values, col_indices, row_ptrs, num_cols)
    }

    /// Returns the column index for a monomial.
    pub fn column_of(&self, m: &PackedMonomial) -> Option<usize> {
        let key = MonomialKey::from(m);
        self.monomial_to_col.get(&key).copied()
    }

    /// Returns the monomial for a column.
    pub fn monomial_of(&self, col: usize) -> Option<&PackedMonomial> {
        self.columns.get(col)
    }

    /// Returns the number of rows.
    pub fn num_rows(&self) -> usize {
        self.rows.len()
    }

    /// Returns the number of columns.
    pub fn num_cols(&self) -> usize {
        self.columns.len()
    }
}

/// Symbolic preprocessing for F4/M5GB.
///
/// Given a set of S-polynomials to reduce, computes the required polynomial
/// multiples to include in the Macaulay matrix.
pub fn symbolic_preprocessing<R: Ring + Clone>(
    spolys: &[(usize, usize, PackedMonomial)],  // (i, j, lcm)
    basis: &[LabeledPoly<R>],
) -> Vec<Vec<PackedMonomial>> {
    let mut multiples: Vec<Vec<PackedMonomial>> = vec![Vec::new(); basis.len()];
    let mut needed: FxHashMap<MonomialKey, ()> = FxHashMap::default();

    // Initial monomials from S-polynomials
    for &(i, j, ref lcm) in spolys {
        if let Some(lm_i) = basis[i].leading_monomial() {
            if let Some(mult) = lcm.div(lm_i) {
                multiples[i].push(mult);
            }
        }
        if let Some(lm_j) = basis[j].leading_monomial() {
            if let Some(mult) = lcm.div(lm_j) {
                multiples[j].push(mult);
            }
        }
    }

    // Collect all monomials that will appear in S-polynomials
    for (poly_idx, mults) in multiples.iter().enumerate() {
        for mult in mults {
            for (_, m) in basis[poly_idx].terms() {
                let product = m.mul(mult);
                let key = MonomialKey::from(&product);
                needed.insert(key, ());
            }
        }
    }

    // Add reducers for each monomial
    let mut processed: FxHashMap<MonomialKey, ()> = FxHashMap::default();

    loop {
        let mut new_needed: Vec<(usize, PackedMonomial)> = Vec::new();

        for (key, _) in needed.iter() {
            if processed.contains_key(key) {
                continue;
            }
            processed.insert(*key, ());

            // Reconstruct monomial from key (simplified)
            let mono = PackedMonomial::new(&key.0[..]);

            // Find a reducer
            for (poly_idx, poly) in basis.iter().enumerate() {
                if let Some(lm) = poly.leading_monomial() {
                    if mono.is_divisible_by(lm) {
                        if let Some(mult) = mono.div(lm) {
                            if !multiples[poly_idx].contains(&mult) {
                                new_needed.push((poly_idx, mult));
                            }
                        }
                        break;
                    }
                }
            }
        }

        if new_needed.is_empty() {
            break;
        }

        // Add new multiples and their monomials
        for (poly_idx, mult) in new_needed {
            multiples[poly_idx].push(mult.clone());

            for (_, m) in basis[poly_idx].terms() {
                let product = m.mul(&mult);
                let key = MonomialKey::from(&product);
                needed.insert(key, ());
            }
        }
    }

    multiples
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signature::Signature;
    use tertius_rings::finite_field::FiniteField;

    type GF101 = FiniteField<101>;

    fn ff(n: u64) -> GF101 {
        GF101::new(n % 101)
    }

    fn make_poly(terms: &[(u64, &[u16])]) -> LabeledPoly<GF101> {
        let terms: Vec<_> = terms
            .iter()
            .map(|(c, exps)| (ff(*c), PackedMonomial::new(exps)))
            .collect();
        let deg = terms.first().map(|(_, m)| m.total_degree()).unwrap_or(0);
        LabeledPoly::new(terms, Signature::generator(0, 3), deg)
    }

    #[test]
    fn test_macaulay_construction() {
        // f = x^2 + y
        let f = make_poly(&[(1, &[2, 0]), (1, &[0, 1])]);
        // g = xy + 1
        let g = make_poly(&[(1, &[1, 1]), (1, &[0, 0])]);

        let polys = vec![f, g];
        let multiples = vec![
            vec![PackedMonomial::one(2)], // f * 1
            vec![PackedMonomial::one(2)], // g * 1
        ];

        let matrix = MacaulayMatrix::construct(&polys, &multiples);

        // Should have 2 rows
        assert_eq!(matrix.num_rows(), 2);

        // Columns should be: x^2, xy, y, 1 (sorted by grevlex)
        assert!(matrix.num_cols() >= 3);
    }
}
