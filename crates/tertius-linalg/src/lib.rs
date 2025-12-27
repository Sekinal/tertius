//! # tertius-linalg
//!
//! High-performance sparse linear algebra for Tertius CAS.
//!
//! This crate provides:
//! - Sparse matrices in CSR (Compressed Sparse Row) format
//! - Dense matrices for small sizes
//! - Parallel matrix operations via rayon
//! - Block Wiedemann algorithm for massive sparse systems
//! - Smith Normal Form computation
//! - SIMD-accelerated GF(2) operations
//!
//! ## Algorithm Selection
//!
//! The crate automatically selects optimal algorithms:
//! - Small matrices (< 100 rows): Dense Gaussian elimination
//! - Large sparse matrices: Block Wiedemann iterative solver
//! - GF(2) matrices: Bit-sliced SIMD operations

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod dense_matrix;
pub mod parallel;
pub mod simd;
pub mod sparse_matrix;

mod berlekamp_massey;
mod block_wiedemann;
mod smith_normal_form;

pub use berlekamp_massey::berlekamp_massey;
pub use block_wiedemann::{BlockWiedemann, BlockWiedemannConfig};
pub use dense_matrix::DenseMatrix;
pub use simd::{Gf2x64, GfPx64};
pub use smith_normal_form::{SmithNormalForm, smith_normal_form};
pub use sparse_matrix::CsrMatrix;

#[cfg(test)]
mod tests;
