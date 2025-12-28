//! # tertius-rings
//!
//! Algebraic structures for Tertius CAS.
//!
//! This crate provides:
//! - Abstract traits: `Ring`, `EuclideanDomain`, `Field`
//! - Concrete implementations: Z, Q, Z_p
//! - Algebraic number fields Q(α)
//! - Polynomial rings R[x]
//!
//! ## Trait Hierarchy
//!
//! ```text
//! Ring
//!  ├── EuclideanDomain
//!  │    └── Field
//!  └── (other specializations)
//! ```

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]

pub mod algebraic;
pub mod finite_field;
pub mod integers;
pub mod poly_ring;
pub mod rationals;
pub mod splitting_field;
pub mod traits;

pub use algebraic::{AlgebraicField, AlgebraicNumber};
pub use finite_field::FiniteField;
pub use integers::Z;
pub use rationals::Q;
pub use splitting_field::SplittingField;
pub use traits::{EuclideanDomain, Field, Ring};
