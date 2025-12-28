//! Power series expansions for the Tertius computer algebra system.
//!
//! This crate provides:
//! - [`PowerSeries`]: Lazy evaluation with memoization
//! - Taylor expansions around a point
//! - Laurent series (with negative powers)
//! - Puiseux series (with fractional powers)
//!
//! # Key algorithms
//!
//! - Composition: Brent-Kung O(√n M(n))
//! - Inversion: Newton iteration
//! - Exp/Log: Fast power series algorithms

pub mod power_series;
pub mod taylor;
pub mod laurent;
pub mod ops;

pub use power_series::{PowerSeries, SeriesCoeff};
pub use taylor::TaylorExpansion;
pub use laurent::LaurentSeries;
