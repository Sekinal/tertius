//! # tertius-simplify
//!
//! Equality saturation-based simplification for Tertius CAS.
//!
//! This crate uses the `egg` library to provide:
//! - E-graph based term rewriting
//! - Algebraic simplification rules
//! - Trigonometric identities
//! - Custom cost functions for extraction
//!
//! ## Equality Saturation vs. Greedy Rewriting
//!
//! Unlike greedy term rewriting, equality saturation explores all
//! possible rewrite paths simultaneously, avoiding local minima
//! and the phase ordering problem.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]

pub mod cost;
pub mod engine;
pub mod language;
pub mod rules;

pub use engine::Simplifier;
pub use language::TertiusLang;
