//! # Tertius
//!
//! A next-generation Computer Algebra System written in Rust.
//!
//! Tertius aims to provide Mathematica-level symbolic computation
//! with the performance and safety guarantees of Rust.
//!
//! ## Features
//!
//! - **High-Performance Core**: Arena-allocated DAG with hash-consing
//! - **Arbitrary Precision**: Full support for big integers and rationals
//! - **Algebraic Structures**: Rings, fields, algebraic extensions
//! - **Polynomial Arithmetic**: Dense/sparse with Karatsuba, FFT, Geobuckets
//! - **Smart Simplification**: Equality saturation via e-graphs
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use tertius::prelude::*;
//!
//! let mut ctx = Context::new();
//! let x = ctx.symbol("x");
//! let expr = (x + 1).pow(2);
//! let expanded = ctx.expand(expr);
//! ```

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]

pub use tertius_core as core;
pub use tertius_integers as integers;
pub use tertius_poly as poly;
pub use tertius_rings as rings;
pub use tertius_simplify as simplify;

/// Prelude module for convenient imports.
pub mod prelude {
    pub use tertius_core::{ExprArena, ExprHandle, ExprNode};
    pub use tertius_integers::{Integer, Rational};
    pub use tertius_poly::{DensePoly, SparsePoly};
    pub use tertius_rings::{Field, Ring, Q, Z};
    pub use tertius_simplify::Simplifier;
}
