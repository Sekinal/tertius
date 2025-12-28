//! Symbolic limit computation for the Tertius computer algebra system.
//!
//! This crate implements the Gruntz algorithm for computing limits of
//! expressions as a variable approaches infinity (or other values).
//!
//! # Algorithm Overview
//!
//! The Gruntz algorithm works by:
//! 1. Finding the Most Rapidly Varying (MRV) subexpressions
//! 2. Rewriting the expression in terms of a single dominant subexpression
//! 3. Computing a series expansion
//! 4. Extracting the limit from the leading term
//!
//! # Example
//!
//! ```ignore
//! use tertius_limits::limit;
//!
//! // Compute lim (x→∞) (x + 1)/x = 1
//! let result = limit("(x + 1)/x", "x", Infinity);
//! ```

pub mod comparison;
pub mod gruntz;
pub mod mrv;

pub use comparison::ComparisonClass;
pub use gruntz::{compute_limit, Limit, LimitResult};
pub use mrv::{MrvSet, find_mrv};
