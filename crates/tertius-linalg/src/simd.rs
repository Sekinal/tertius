//! SIMD-accelerated operations for finite fields.
//!
//! This module provides bit-sliced representations for efficient
//! parallel computation over GF(2) and small prime fields.
//!
//! # GF(2) Bit-Slicing
//!
//! For GF(2), we pack 64 field elements into a single u64.
//! Addition is XOR, multiplication is AND. This allows processing
//! 64 field operations with a single CPU instruction.
//!
//! # Small Prime Fields
//!
//! For small primes that fit in u8 or u16, we can pack multiple
//! elements into wider registers for SIMD processing.

use std::ops::{Add, BitAnd, BitXor, Mul, Not};

use num_traits::{One, Zero};
use rayon::prelude::*;

/// 64 GF(2) elements packed into a u64.
///
/// Each bit represents one field element:
/// - 0 = zero
/// - 1 = one
///
/// Addition (XOR) and multiplication (AND) operate on all 64 elements
/// simultaneously.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct Gf2x64(pub u64);

impl Gf2x64 {
    /// All zeros.
    pub const ZERO: Self = Self(0);

    /// All ones.
    pub const ONE: Self = Self(u64::MAX);

    /// Creates a new packed value.
    #[inline]
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Creates a value with a single bit set.
    #[inline]
    #[must_use]
    pub const fn single(index: usize) -> Self {
        Self(1u64 << index)
    }

    /// Returns the raw u64 value.
    #[inline]
    #[must_use]
    pub const fn raw(self) -> u64 {
        self.0
    }

    /// Gets the bit at the specified index.
    #[inline]
    #[must_use]
    pub const fn get(self, index: usize) -> bool {
        (self.0 >> index) & 1 == 1
    }

    /// Sets the bit at the specified index.
    #[inline]
    #[must_use]
    pub const fn set(self, index: usize, value: bool) -> Self {
        if value {
            Self(self.0 | (1u64 << index))
        } else {
            Self(self.0 & !(1u64 << index))
        }
    }

    /// Counts the number of set bits (population count).
    #[inline]
    #[must_use]
    pub const fn popcount(self) -> u32 {
        self.0.count_ones()
    }

    /// Returns true if any bit is set.
    #[inline]
    #[must_use]
    pub const fn any(self) -> bool {
        self.0 != 0
    }

    /// Returns true if all bits are zero.
    #[inline]
    #[must_use]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// Computes the dot product of two packed vectors.
    ///
    /// This is equivalent to XOR of (a AND b), then reducing to a single bit.
    #[inline]
    #[must_use]
    pub const fn dot(self, other: Self) -> bool {
        (self.0 & other.0).count_ones() % 2 == 1
    }
}

impl Zero for Gf2x64 {
    fn zero() -> Self {
        Self::ZERO
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl One for Gf2x64 {
    fn one() -> Self {
        Self::ONE
    }
}

impl Add for Gf2x64 {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl Mul for Gf2x64 {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }
}

impl BitXor for Gf2x64 {
    type Output = Self;

    #[inline]
    fn bitxor(self, other: Self) -> Self {
        Self(self.0 ^ other.0)
    }
}

impl BitAnd for Gf2x64 {
    type Output = Self;

    #[inline]
    fn bitand(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }
}

impl Not for Gf2x64 {
    type Output = Self;

    #[inline]
    fn not(self) -> Self {
        Self(!self.0)
    }
}

/// 64 elements of GF(p) for small primes, packed into u64.
///
/// Each element uses `bits_per_element` bits.
/// For p=3, we use 2 bits (values 0, 1, 2).
/// For p=7, we use 3 bits (values 0-6).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GfPx64<const P: u64> {
    /// Packed elements.
    data: u64,
    /// Number of elements stored.
    count: usize,
}

impl<const P: u64> GfPx64<P> {
    /// Bits needed to represent an element of GF(P).
    const BITS_PER_ELEMENT: usize = (64 - P.leading_zeros()) as usize;

    /// Maximum number of elements that fit in a u64.
    const MAX_COUNT: usize = 64 / Self::BITS_PER_ELEMENT;

    /// Mask for a single element.
    const ELEMENT_MASK: u64 = (1u64 << Self::BITS_PER_ELEMENT) - 1;

    /// Creates a new packed value with all zeros.
    #[must_use]
    pub const fn zero() -> Self {
        Self { data: 0, count: 0 }
    }

    /// Creates a packed value from an iterator of elements.
    #[must_use]
    pub fn from_iter(iter: impl IntoIterator<Item = u64>) -> Self {
        let mut data = 0u64;
        let mut count = 0;
        for val in iter {
            if count >= Self::MAX_COUNT {
                break;
            }
            let shift = count * Self::BITS_PER_ELEMENT;
            data |= (val % P) << shift;
            count += 1;
        }
        Self { data, count }
    }

    /// Gets the element at the specified index.
    #[must_use]
    pub fn get(&self, index: usize) -> u64 {
        let shift = index * Self::BITS_PER_ELEMENT;
        (self.data >> shift) & Self::ELEMENT_MASK
    }

    /// Sets the element at the specified index.
    pub fn set(&mut self, index: usize, value: u64) {
        let shift = index * Self::BITS_PER_ELEMENT;
        let mask = Self::ELEMENT_MASK << shift;
        self.data = (self.data & !mask) | ((value % P) << shift);
        if index >= self.count {
            self.count = index + 1;
        }
    }

    /// Returns the number of elements.
    #[must_use]
    pub fn len(&self) -> usize {
        self.count
    }

    /// Returns true if empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Extracts all elements to a vector.
    #[must_use]
    pub fn to_vec(&self) -> Vec<u64> {
        (0..self.count).map(|i| self.get(i)).collect()
    }
}

impl<const P: u64> Add for GfPx64<P> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        // Element-wise addition mod P
        // This is tricky because of potential overflow between packed elements
        // For now, use scalar approach
        let count = self.count.max(other.count);
        let mut result = Self::zero();
        result.count = count;

        for i in 0..count {
            let a = self.get(i);
            let b = other.get(i);
            result.set(i, (a + b) % P);
        }
        result
    }
}

impl<const P: u64> Mul for GfPx64<P> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let count = self.count.max(other.count);
        let mut result = Self::zero();
        result.count = count;

        for i in 0..count {
            let a = self.get(i);
            let b = other.get(i);
            result.set(i, (a * b) % P);
        }
        result
    }
}

/// Bit-sliced matrix over GF(2).
///
/// Each column is stored as a Vec<Gf2x64>, allowing 64 rows
/// to be processed in parallel.
#[derive(Clone, Debug)]
pub struct Gf2Matrix {
    /// Columns of the matrix, each storing 64 rows.
    columns: Vec<Vec<Gf2x64>>,
    /// Number of rows.
    num_rows: usize,
    /// Number of columns.
    num_cols: usize,
}

impl Gf2Matrix {
    /// Creates a new zero matrix.
    #[must_use]
    pub fn zeros(num_rows: usize, num_cols: usize) -> Self {
        let num_chunks = (num_rows + 63) / 64;
        Self {
            columns: vec![vec![Gf2x64::ZERO; num_chunks]; num_cols],
            num_rows,
            num_cols,
        }
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

    /// Gets the entry at (row, col).
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> bool {
        let chunk = row / 64;
        let bit = row % 64;
        self.columns[col][chunk].get(bit)
    }

    /// Sets the entry at (row, col).
    pub fn set(&mut self, row: usize, col: usize, value: bool) {
        let chunk = row / 64;
        let bit = row % 64;
        self.columns[col][chunk] = self.columns[col][chunk].set(bit, value);
    }

    /// Adds column src to column dst (XOR).
    pub fn add_column(&mut self, dst: usize, src: usize) {
        let num_chunks = (self.num_rows + 63) / 64;
        for chunk in 0..num_chunks {
            self.columns[dst][chunk] = self.columns[dst][chunk] ^ self.columns[src][chunk];
        }
    }

    /// Matrix-vector multiply over GF(2).
    ///
    /// The vector is represented as bits in Gf2x64 values.
    #[must_use]
    pub fn mv(&self, x: &[Gf2x64]) -> Vec<Gf2x64> {
        let num_chunks = (self.num_rows + 63) / 64;
        let mut result = vec![Gf2x64::ZERO; num_chunks];

        for (col_idx, col) in self.columns.iter().enumerate() {
            let x_chunk = col_idx / 64;
            let x_bit = col_idx % 64;

            if x_chunk < x.len() && x[x_chunk].get(x_bit) {
                for (chunk_idx, &col_chunk) in col.iter().enumerate() {
                    result[chunk_idx] = result[chunk_idx] ^ col_chunk;
                }
            }
        }

        result
    }

    /// Parallel matrix-vector multiply over GF(2).
    #[must_use]
    pub fn mv_parallel(&self, x: &[Gf2x64]) -> Vec<Gf2x64> {
        let num_chunks = (self.num_rows + 63) / 64;

        (0..num_chunks)
            .into_par_iter()
            .map(|chunk_idx| {
                let mut result = Gf2x64::ZERO;
                for (col_idx, col) in self.columns.iter().enumerate() {
                    let x_chunk = col_idx / 64;
                    let x_bit = col_idx % 64;

                    if x_chunk < x.len() && x[x_chunk].get(x_bit) {
                        result = result ^ col[chunk_idx];
                    }
                }
                result
            })
            .collect()
    }

    /// Row reduction (Gaussian elimination) over GF(2).
    ///
    /// Modifies the matrix in place and returns the rank.
    pub fn row_reduce(&mut self) -> usize {
        let mut pivot_row = 0;

        for col in 0..self.num_cols {
            // Find pivot
            let mut found = false;
            for row in pivot_row..self.num_rows {
                if self.get(row, col) {
                    // Swap rows if needed
                    if row != pivot_row {
                        self.swap_rows(pivot_row, row);
                    }
                    found = true;
                    break;
                }
            }

            if !found {
                continue;
            }

            // Eliminate other entries in this column
            for row in 0..self.num_rows {
                if row != pivot_row && self.get(row, col) {
                    self.add_row(row, pivot_row);
                }
            }

            pivot_row += 1;
        }

        pivot_row // This is the rank
    }

    /// Swaps two rows.
    fn swap_rows(&mut self, i: usize, j: usize) {
        let chunk_i = i / 64;
        let bit_i = i % 64;
        let chunk_j = j / 64;
        let bit_j = j % 64;

        for col in &mut self.columns {
            let val_i = col[chunk_i].get(bit_i);
            let val_j = col[chunk_j].get(bit_j);
            col[chunk_i] = col[chunk_i].set(bit_i, val_j);
            col[chunk_j] = col[chunk_j].set(bit_j, val_i);
        }
    }

    /// Adds row src to row dst.
    fn add_row(&mut self, dst: usize, src: usize) {
        let chunk_dst = dst / 64;
        let bit_dst = dst % 64;
        let chunk_src = src / 64;
        let bit_src = src % 64;

        for col in &mut self.columns {
            let val_src = col[chunk_src].get(bit_src);
            let val_dst = col[chunk_dst].get(bit_dst);
            col[chunk_dst] = col[chunk_dst].set(bit_dst, val_dst ^ val_src);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gf2x64_basic() {
        let a = Gf2x64::new(0b1010);
        let b = Gf2x64::new(0b1100);

        // Addition (XOR)
        assert_eq!((a + b).raw(), 0b0110);

        // Multiplication (AND)
        assert_eq!((a * b).raw(), 0b1000);
    }

    #[test]
    fn test_gf2x64_get_set() {
        let mut v = Gf2x64::ZERO;
        assert!(!v.get(5));

        v = v.set(5, true);
        assert!(v.get(5));
        assert!(!v.get(4));

        v = v.set(5, false);
        assert!(!v.get(5));
    }

    #[test]
    fn test_gf2x64_popcount() {
        assert_eq!(Gf2x64::new(0b1010101).popcount(), 4);
        assert_eq!(Gf2x64::ZERO.popcount(), 0);
        assert_eq!(Gf2x64::ONE.popcount(), 64);
    }

    #[test]
    fn test_gf2x64_dot() {
        let a = Gf2x64::new(0b1010);
        let b = Gf2x64::new(0b1100);
        // AND: 0b1000, popcount = 1, so dot = true (odd)
        assert!(a.dot(b));

        let c = Gf2x64::new(0b1111);
        let d = Gf2x64::new(0b1111);
        // AND: 0b1111, popcount = 4, so dot = false (even)
        assert!(!c.dot(d));
    }

    #[test]
    fn test_gfpx64_basic() {
        let mut v: GfPx64<7> = GfPx64::zero();
        v.set(0, 3);
        v.set(1, 5);
        v.set(2, 10); // 10 % 7 = 3

        assert_eq!(v.get(0), 3);
        assert_eq!(v.get(1), 5);
        assert_eq!(v.get(2), 3);
    }

    #[test]
    fn test_gf2_matrix() {
        let mut m = Gf2Matrix::zeros(4, 4);

        // Set identity
        for i in 0..4 {
            m.set(i, i, true);
        }

        // Check diagonal
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m.get(i, j), i == j);
            }
        }
    }

    #[test]
    fn test_gf2_matrix_mv() {
        // 2x2 matrix [[1,1],[0,1]]
        let mut m = Gf2Matrix::zeros(2, 2);
        m.set(0, 0, true);
        m.set(0, 1, true);
        m.set(1, 1, true);

        // Vector [1, 0]
        let x = vec![Gf2x64::new(0b01)];
        let y = m.mv(&x);

        // Result should be [1, 0]
        assert!(y[0].get(0));
        assert!(!y[0].get(1));
    }

    #[test]
    fn test_gf2_matrix_row_reduce() {
        // Matrix [[1,0,1],[1,1,0],[0,1,1]] - has rank 2 in GF(2)
        // det = 1*(1*1+0*1) + 0 + 1*(1*1+1*0) = 1 + 1 = 0 in GF(2)
        let mut m = Gf2Matrix::zeros(3, 3);
        m.set(0, 0, true);
        m.set(0, 2, true);
        m.set(1, 0, true);
        m.set(1, 1, true);
        m.set(2, 1, true);
        m.set(2, 2, true);

        let rank = m.row_reduce();
        assert_eq!(rank, 2); // Rank 2 (singular in GF(2))

        // Test a full-rank matrix: [[1,0,0],[0,1,0],[0,0,1]] (identity)
        let mut id = Gf2Matrix::zeros(3, 3);
        id.set(0, 0, true);
        id.set(1, 1, true);
        id.set(2, 2, true);

        let rank_id = id.row_reduce();
        assert_eq!(rank_id, 3); // Full rank
    }
}
