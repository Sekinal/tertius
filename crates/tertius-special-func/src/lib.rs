//! Special Functions for Tertius CAS
//!
//! This crate provides representations and transformations for special functions
//! that arise in symbolic integration:
//!
//! - **Polylogarithms**: Li_n(x) = Σ x^k/k^n
//! - **Error function**: erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
//! - **Logarithmic integral**: li(x) = ∫₀ˣ dt/ln(t)
//! - **Elliptic integrals**: F(φ,k), E(φ,k), Π(n;φ,k)
//! - **Hypergeometric functions**: ₂F₁(a,b;c;z)
//!
//! # Purpose
//!
//! When the Risch algorithm determines an integral is non-elementary,
//! these special functions may provide an answer. For example:
//!
//! - ∫ e^(-x²) dx = (√π/2) erf(x)
//! - ∫ 1/ln(x) dx = li(x)
//! - ∫ ln(x)/(x-1) dx = Li₂(x) + C

pub mod polylog;
pub mod error_func;
pub mod elliptic;
pub mod hypergeometric;
pub mod catalog;

pub use polylog::{Polylogarithm, dilog, trilog};
pub use error_func::{ErrorFunction, error_func_series};
pub use elliptic::{EllipticF, EllipticE, EllipticPi};
pub use hypergeometric::{Hypergeometric2F1, hypergeometric_series};
pub use catalog::{SpecialFunction, recognize_special_function};
