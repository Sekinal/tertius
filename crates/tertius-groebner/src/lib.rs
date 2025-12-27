//! M5GB Gröbner basis algorithm.
//!
//! This crate provides a hybrid algorithm combining:
//! - F5 signatures for early detection of useless pairs
//! - M4GB-style reductor caching for efficient reductions
//! - Parallel matrix reduction via rayon

pub mod criteria;
pub mod labeled_poly;
pub mod m5gb;
pub mod macaulay;
pub mod monomial;
pub mod parallel_reduce;
pub mod reductor_store;
pub mod signature;

pub use m5gb::M5GB;
