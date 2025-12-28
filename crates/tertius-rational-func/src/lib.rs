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

mod rational_func;
mod arithmetic;
pub mod partial_fractions;
pub mod hermite;

pub use rational_func::RationalFunction;
