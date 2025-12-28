//! Polynomial factorization algorithms for Tertius CAS.
//!
//! This crate provides:
//! - **LLL**: Lattice reduction for Van Hoeij's algorithm
//! - **Berlekamp-Zassenhaus**: Factorization over finite fields
//! - **Cantor-Zassenhaus**: Probabilistic equal-degree factorization
//! - **Hensel Lifting**: Lifting factors from Z_p to Z
//! - **Van Hoeij**: Lattice-based univariate factorization over Z
//! - **Lecerf**: Multivariate factorization
//!
//! # Parallelism
//!
//! All algorithms use rayon for parallel execution where applicable.

pub mod berlekamp;
pub mod cantor_zassenhaus;
pub mod hensel;
pub mod leading_coeff;
pub mod lecerf;
pub mod lll;
pub mod multivariate;
pub mod multivariate_hensel;
pub mod squarefree;
pub mod univariate;

// Re-exports
pub use berlekamp::berlekamp_factor;
pub use cantor_zassenhaus::cantor_zassenhaus_factor;
pub use hensel::hensel_lift;
pub use leading_coeff::{compute_leading_coeff_info, LeadingCoeffInfo};
pub use lecerf::{lecerf_factor_multivariate, LecerfMultivariateResult};
pub use lll::lll_reduce;
pub use multivariate::{
    factor_bivariate, factor_multivariate, factor_multivariate_batch, lecerf_factor,
    MultivariateFactorResult,
};
pub use multivariate_hensel::{multivariate_hensel_lift, MultivariateHenselResult};
pub use squarefree::squarefree_factorization;
pub use univariate::van_hoeij_factor;
