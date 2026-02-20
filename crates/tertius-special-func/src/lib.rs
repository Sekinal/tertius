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

pub mod catalog;
pub mod elliptic;
pub mod error_func;
pub mod hypergeometric;
pub mod physics;
pub mod polylog;

pub use catalog::{recognize_special_function, SpecialFunction};
pub use elliptic::{EllipticE, EllipticF, EllipticPi};
pub use error_func::{error_func_series, ErrorFunction};
pub use hypergeometric::{hypergeometric_series, Hypergeometric2F1};
pub use physics::{
    bessel_j, derivative_rule, hermite_h, integral_rule, laguerre_l, legendre_p, simplify_rule,
    spherical_harmonic_y_00, PhysicsSpecialExpr,
};
pub use polylog::{dilog, trilog, Polylogarithm};
