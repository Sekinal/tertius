//! # tertius-poly
//!
//! High-performance polynomial arithmetic for Tertius CAS.
//!
//! This crate provides:
//! - Dense univariate polynomials with Karatsuba/FFT multiplication
//! - Sparse multivariate polynomials with bit-packed monomials
//! - Geobucket-based sparse multiplication
//! - Polynomial GCD algorithms
//!
//! ## Algorithm Selection
//!
//! Multiplication automatically selects the optimal algorithm:
//! - Degree < 32: Schoolbook O(n²)
//! - Degree 32-1023: Karatsuba O(n^1.58)
//! - Degree >= 1024: FFT/NTT O(n log n)

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]

pub mod algorithms;
pub mod dense;
pub mod monomial;
pub mod ordering;
pub mod sparse;

#[cfg(test)]
mod proptests;

pub use dense::DensePoly;
pub use monomial::PackedMonomial;
pub use ordering::MonomialOrder;
pub use sparse::SparsePoly;
