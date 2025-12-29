//! # tertius-numtheory
//!
//! Number theory algorithms for Tertius CAS.
//!
//! This crate provides:
//! - **Primality testing**: Trial division, Miller-Rabin, Baillie-PSW (deterministic for 64-bit)
//! - **Integer factorization**: Trial division, Pollard's rho, p-1, p+1, ECM, SIQS
//! - **Arithmetic functions**: Euler's totient φ(n), divisor functions σ_k(n), Möbius μ(n)
//! - **Modular arithmetic**: CRT, Legendre/Jacobi symbols, Tonelli-Shanks square roots
//! - **Prime generation**: Segmented sieve, prime counting π(x)
//!
//! ## Algorithm Selection
//!
//! The factorization dispatcher automatically selects the best algorithm based on input size:
//!
//! | Size | Algorithm | Performance |
//! |------|-----------|-------------|
//! | < 10^6 | Trial division | Instant |
//! | 10^6 - 10^15 | Pollard's rho | < 1 second |
//! | 10^15 - 10^25 | ECM | Seconds to minutes |
//! | 10^25 - 10^100 | SIQS | Minutes to hours |
//!
//! ## Example
//!
//! ```
//! use tertius_numtheory::{factor, is_prime, euler_phi};
//! use tertius_integers::Integer;
//!
//! // Primality testing
//! assert!(is_prime(&Integer::from(1000000007i64)));
//!
//! // Factorization
//! let n = Integer::from(1234567890i64);
//! let factors = factor(&n);
//! // factors = [(2, 1), (3, 2), (5, 1), (3607, 1), (3803, 1)]
//!
//! // Euler's totient
//! let phi = euler_phi(&Integer::from(12i64));
//! assert_eq!(phi, Integer::from(4i64)); // φ(12) = 4
//! ```

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::module_name_repetitions)]

pub mod primality;
pub mod factorization;
pub mod arithmetic;
pub mod congruences;
pub mod primes;

// Re-export main types and functions
pub use primality::{is_prime, is_probable_prime, PrimalityResult};
pub use factorization::{factor, Factorization};
pub use arithmetic::{euler_phi, carmichael_lambda, divisors, divisor_count, divisor_sum, mobius, radical};
pub use congruences::{crt, legendre_symbol, jacobi_symbol, tonelli_shanks};
pub use primes::{sieve_of_eratosthenes, nth_prime, prime_count, primes_up_to};
