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
pub mod special_forms;
pub mod numerical;
pub mod definite;

// Re-exports
pub use rational::{
    integrate_rational, integrate_rational_q, integrate_rational_with_algebraic, AlgebraicIntegrationResult,
    RationalIntegrationResult,
};
pub use rothstein_trager::{AlgebraicLogarithmicPart, LogarithmicPart, rothstein_trager_q};
pub use risch::{risch_integrate, IntegralExpression, RischResult};
pub use proofs::{prove_non_elementary, NonIntegrabilityProof, NonIntegrabilityReason};
pub use special_forms::{recognize_special_form, SpecialFormResult, SpecialFunction, SpecialIntegral};
