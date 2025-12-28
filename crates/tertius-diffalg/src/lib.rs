//! Differential Algebra Framework for the Risch Algorithm
//!
//! This crate provides the foundational structures for symbolic integration:
//!
//! - **Derivations**: The abstract differentiation operator D
//! - **Transcendental Tower**: K ⊂ K(θ₁) ⊂ K(θ₂) ⊂ ... where each θᵢ is
//!   logarithmic or exponential
//! - **Tower Elements**: Representations of elements in the tower
//!
//! # The Risch Algorithm Framework
//!
//! The Risch algorithm works by building a tower of transcendental extensions:
//!
//! 1. Start with base field K = Q(x) (rational functions)
//! 2. Add transcendentals θ = log(u) or θ = exp(u)
//! 3. Elements are rational functions in θ with coefficients in previous level
//! 4. Recursively integrate down the tower
//!
//! # Key Properties
//!
//! For θ = log(u): θ' = u'/u
//! For θ = exp(u): θ' = u'·θ

pub mod derivation;
pub mod tower;
pub mod element;

pub use derivation::Derivation;
pub use tower::{TranscendentalTower, TranscendentalType, TowerLevel};
pub use element::TowerElement;
