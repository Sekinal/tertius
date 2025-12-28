//! Non-integrability proofs for the Risch algorithm.
//!
//! This module provides formal proofs that certain integrals have no
//! elementary antiderivative, based on Liouville's theorem and the
//! structure theory of the Risch algorithm.
//!
//! # Liouville's Theorem
//!
//! If f ∈ K(θ) where θ is transcendental over K with θ' = θu'/u (logarithmic)
//! or θ' = θu' (exponential), then:
//!
//! ∫f dx is elementary iff ∫f dx = v₀ + Σ cᵢ log(vᵢ)
//!
//! where v₀, vᵢ ∈ K(θ) and cᵢ are constants.
//!
//! # Proof Structure
//!
//! A non-integrability proof consists of:
//! 1. The integrand and its transcendental tower
//! 2. The obstruction (why the Risch equations have no solution)
//! 3. A human-readable explanation

pub mod liouville;

pub use liouville::{
    NonIntegrabilityProof,
    NonIntegrabilityReason,
    prove_non_elementary,
    check_liouville_obstruction,
};
