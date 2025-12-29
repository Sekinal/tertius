//! Combinatorics for the Tertius CAS.
//!
//! This crate provides combinatorial functions and sequences:
//!
//! - Factorials and binomial coefficients
//! - Stirling numbers (first and second kind)
//! - Bell numbers
//! - Catalan numbers
//! - Bernoulli numbers
//! - Integer partitions
//!
//! # Examples
//!
//! ```
//! use tertius_combinatorics::{factorial, binomial, catalan, bell};
//! use tertius_integers::Integer;
//!
//! // 5! = 120
//! assert_eq!(factorial(5), Integer::from(120i64));
//!
//! // C(10, 3) = 120
//! assert_eq!(binomial(10, 3), Integer::from(120i64));
//!
//! // C_5 = 42
//! assert_eq!(catalan(5), Integer::from(42i64));
//!
//! // B_5 = 52
//! assert_eq!(bell(5), Integer::from(52i64));
//! ```

mod factorial;
mod binomial;
mod stirling;
mod bell;
mod catalan;
mod bernoulli;
mod partitions;

pub use factorial::{factorial, falling_factorial, rising_factorial, double_factorial};
pub use binomial::binomial;
pub use stirling::{stirling1, stirling2};
pub use bell::bell;
pub use catalan::catalan;
pub use bernoulli::bernoulli;
pub use partitions::{partition_count, partitions};
