//! # tertius-core
//!
//! Core expression engine for the Tertius Computer Algebra System.
//!
//! This crate provides:
//! - Arena-allocated expression storage with hash-consing
//! - Type-safe expression handles
//! - O(1) structural equality via interning
//!
//! ## Design Principles
//!
//! - **Data-Oriented Design**: Expressions stored contiguously in arena for cache efficiency
//! - **Hash-Consing**: Every structurally unique expression stored exactly once
//! - **Zero-Cost Handles**: 32-bit indices instead of pointers

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]

pub mod arena;
pub mod expr;
pub mod handle;
pub mod intern;

pub use arena::ExprArena;
pub use expr::{ExprNode, FunctionId, SymbolId};
pub use handle::ExprHandle;
