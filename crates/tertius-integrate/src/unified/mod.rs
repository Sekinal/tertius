//! Unified Integration API
//!
//! This module provides a Mathematica-style `integrate()` function that
//! automatically handles any expression type by classifying the integrand
//! and dispatching to the appropriate backend.
//!
//! # Example
//!
//! ```ignore
//! use tertius_integrate::integrate;
//! use tertius_core::ExprArena;
//!
//! let mut arena = ExprArena::new();
//! let x = arena.symbol("x");
//! let x_squared = arena.pow(x, arena.integer(2));
//!
//! let result = integrate(&mut arena, x_squared, x);
//! // result contains x^3/3
//! ```

pub mod result;
pub mod classifier;
pub mod convert;
pub mod dispatch;
pub mod api;

pub use api::{integrate, integrate_definite, integrate_with_options};
pub use result::{
    IntegrationMethod, IntegrationResult, NumericalResult, SymbolicAntiderivative,
    SpecialFunctionResult, EllipticResult, NonElementaryResult, UnknownReason,
};
pub use dispatch::IntegrationOptions;
