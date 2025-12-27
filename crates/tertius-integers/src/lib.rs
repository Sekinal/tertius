//! # tertius-integers
//!
//! Arbitrary precision integer and rational arithmetic for Tertius CAS.
//!
//! This crate wraps `dashu` to provide:
//! - Arbitrary precision integers (`Integer`)
//! - Arbitrary precision rationals (`Rational`)
//! - Modular arithmetic (`ModInt`)
//!
//! ## Performance Notes
//!
//! - Small integers (fitting in a machine word) use stack allocation
//! - Large integers are heap-allocated with GMP-like performance

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]

pub mod integer;
pub mod modular;
pub mod rational;

#[cfg(test)]
mod proptests;

pub use integer::Integer;
pub use modular::{ModInt, ModIntDyn};
pub use rational::Rational;
