//! Symbolic differentiation for Tertius CAS.
//!
//! This crate provides symbolic differentiation of expressions,
//! including partial derivatives, gradients, Jacobians, and Hessians.
//!
//! # Examples
//!
//! ```
//! use tertius_core::arena::ExprArena;
//! use tertius_diff::{diff, diff_n};
//!
//! let mut arena = ExprArena::new();
//! let x = arena.symbol("x");
//!
//! // Create x^2
//! let two = arena.integer(2);
//! let x_squared = arena.pow(x, two);
//!
//! // Differentiate: d/dx (x^2) = 2x
//! let derivative = diff(&mut arena, x_squared, x);
//! ```

mod diff;
mod partial;
mod rules;

pub use diff::{diff, diff_n};
pub use partial::{gradient, hessian, jacobian, partial};
