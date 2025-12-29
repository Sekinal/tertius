//! Symbolic Integration Engine for Tertius CAS
//!
//! This crate provides a unified integration API like Mathematica's `Integrate[]`,
//! automatically dispatching to appropriate algorithms based on the integrand structure.
//!
//! # Quick Start
//!
//! ```ignore
//! use tertius_integrate::{integrate, integrate_definite, integrate_numerical};
//! use tertius_core::ExprArena;
//!
//! let mut arena = ExprArena::new();
//! let x = arena.symbol("x");
//!
//! // Indefinite: ∫ x² dx = x³/3
//! let x_squared = arena.pow(x, arena.integer(2));
//! let result = integrate(&mut arena, x_squared, x);
//!
//! // Definite: ∫₀¹ x² dx = 1/3
//! let result = integrate_definite(&mut arena, x_squared, x, arena.integer(0), arena.integer(1));
//!
//! // Numerical: ∫₀^π sin(x) dx = 2
//! let result = integrate_numerical(&|x| x.sin(), 0.0, std::f64::consts::PI, 1e-10);
//! ```
//!
//! # Features
//!
//! - **Unified API**: Single `integrate()` function handles any expression
//! - **Automatic dispatch**: Classifies integrand and routes to appropriate algorithm
//! - **Rational functions**: Hermite reduction + Rothstein-Trager
//! - **Transcendental extensions**: Risch algorithm for log/exp towers
//! - **Algebraic functions**: Genus-based classification and elliptic integrals
//! - **Special functions**: Recognition of erf, Ei, Li, etc.
//! - **Numerical fallback**: Adaptive Gauss-Kronrod quadrature
//! - **Non-integrability proofs**: Formal proofs when no elementary antiderivative exists

pub mod rational;
pub mod rothstein_trager;
pub mod polynomial;
pub mod risch;
pub mod proofs;
pub mod special_forms;
pub mod numerical;
pub mod definite;
pub mod algebraic;
pub mod special_definite;
pub mod unified;

// Primary API - unified integration
pub use unified::{
    integrate, integrate_definite, integrate_with_options,
    IntegrationResult, IntegrationMethod, IntegrationOptions,
    SymbolicAntiderivative, NumericalResult, SpecialFunctionResult,
    EllipticResult, NonElementaryResult, UnknownReason,
};
pub use unified::api::{integrate_numerical, integrate_numerical_with_params};

// Re-exports for direct backend access
pub use rational::{
    integrate_rational, integrate_rational_q, integrate_rational_with_algebraic, AlgebraicIntegrationResult,
    RationalIntegrationResult,
};
pub use rothstein_trager::{AlgebraicLogarithmicPart, LogarithmicPart, rothstein_trager_q};
pub use risch::{risch_integrate, IntegralExpression, RischResult};
pub use proofs::{prove_non_elementary, NonIntegrabilityProof, NonIntegrabilityReason as ProofReason};
pub use special_forms::{recognize_special_form, SpecialFormResult, SpecialFunction as SpecialFormFunc, SpecialIntegral};
