//! Polynomial algorithms.
//!
//! This module contains high-performance implementations of:
//! - Karatsuba multiplication
//! - FFT/NTT-based multiplication
//! - Kronecker substitution
//! - Geobucket-based sparse multiplication
//! - Polynomial GCD
//! - Berlekamp-Massey algorithm
//! - Ben-Or/Tiwari sparse interpolation
//! - Sparse GCD (Hu-Monagan)
//! - Resultants

pub mod fft;
pub mod gcd;
pub mod geobucket;
pub mod karatsuba;
pub mod ntt;

pub mod ben_or_tiwari;
pub mod berlekamp_massey;
pub mod resultant;
pub mod sparse_gcd;
