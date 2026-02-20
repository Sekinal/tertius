//! Rational functions P(x)/Q(x) for Tertius CAS.
//!
//! This crate provides:
//! - [`RationalFunction`] type for representing P(x)/Q(x)
//! - Arithmetic operations (add, sub, mul, div)
//! - Partial fraction decomposition
//! - Hermite reduction (reduces to simple poles)
//!
//! These are foundational components for the Risch algorithm
//! for symbolic integration.

mod arithmetic;
pub mod hermite;
pub mod partial_fractions;
mod rational_func;

pub use rational_func::RationalFunction;
