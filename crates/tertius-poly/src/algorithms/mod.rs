//! Polynomial algorithms.
//!
//! This module contains high-performance implementations of:
//! - Karatsuba multiplication
//! - FFT/NTT-based multiplication
//! - Kronecker substitution
//! - Geobucket-based sparse multiplication
//! - Polynomial GCD

pub mod fft;
pub mod gcd;
pub mod karatsuba;
pub mod ntt;

// TODO: Implement these modules
// pub mod kronecker;
// pub mod geobucket;
