//! Symbolic Integration Engine for Tertius CAS
//!
//! This crate implements the Risch algorithm for computing symbolic integrals.
//!
//! # Features
//!
//! - **Rational function integration**: Complete algorithm using Hermite reduction
//!   and Rothstein-Trager for logarithmic parts
//! - **Transcendental extensions**: Integration over log and exp towers
//! - **Special functions**: Recognition and use of polylogarithms, error functions, etc.
//! - **Non-integrability proofs**: Formal proofs when no elementary antiderivative exists
//!
//! # Example
//!
//! ```ignore
//! use tertius_integrate::{integrate_rational, IntegrationResult};
//! use tertius_rational_func::RationalFunction;
//!
//! // Integrate 1/(x^2 - 1)
//! let rf = RationalFunction::new(...);
//! let result = integrate_rational(&rf);
//! ```

pub mod rational;
pub mod rothstein_trager;
pub mod polynomial;
pub mod risch;
pub mod proofs;

// Re-exports
pub use rational::{integrate_rational, RationalIntegrationResult};
pub use rothstein_trager::LogarithmicPart;
pub use risch::{RischResult, IntegralExpression, risch_integrate};
pub use proofs::{NonIntegrabilityProof, NonIntegrabilityReason, prove_non_elementary};
