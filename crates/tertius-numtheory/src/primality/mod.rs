//! Primality testing algorithms.
//!
//! This module provides several primality tests with different trade-offs:
//!
//! - **Trial division**: Exact for small numbers, O(√n)
//! - **Miller-Rabin**: Probabilistic, O(k log³n) where k is the number of witnesses
//! - **Baillie-PSW**: Deterministic for n < 2^64, no known pseudoprimes above that
//!
//! The main entry point is [`is_prime`] which uses Baillie-PSW.

mod trial_division;
mod miller_rabin;
mod lucas;
mod baillie_psw;

pub use trial_division::trial_division;
pub use miller_rabin::{miller_rabin, miller_rabin_deterministic};
pub use lucas::{lucas_test, strong_lucas_test};
pub use baillie_psw::baillie_psw;

use tertius_integers::Integer;

/// Result of a primality test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimalityResult {
    /// The number is definitely prime.
    Prime,
    /// The number is definitely composite.
    Composite,
    /// The number is probably prime (probabilistic test).
    ProbablyPrime,
}

impl PrimalityResult {
    /// Returns true if the result indicates the number is (probably) prime.
    pub fn is_prime(self) -> bool {
        matches!(self, PrimalityResult::Prime | PrimalityResult::ProbablyPrime)
    }
}

/// Tests if n is prime using the Baillie-PSW test.
///
/// This is deterministic for n < 2^64 (no known pseudoprimes).
/// For larger numbers, it is probabilistic but extremely reliable.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::is_prime;
/// use tertius_integers::Integer;
///
/// assert!(is_prime(&Integer::from(2i64)));
/// assert!(is_prime(&Integer::from(1000000007i64)));
/// assert!(!is_prime(&Integer::from(1000000008i64)));
/// ```
pub fn is_prime(n: &Integer) -> bool {
    baillie_psw(n).is_prime()
}

/// Tests if n is probably prime with error probability < 4^(-k).
///
/// Uses k rounds of Miller-Rabin with random witnesses.
/// For most applications, k=25 is sufficient.
pub fn is_probable_prime(n: &Integer, k: u32) -> PrimalityResult {
    // First do cheap trial division
    if let Some(result) = trial_division::quick_check(n) {
        return result;
    }

    // Then do k rounds of Miller-Rabin
    if miller_rabin(n, k) {
        PrimalityResult::ProbablyPrime
    } else {
        PrimalityResult::Composite
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_primes() {
        let small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        for &p in &small_primes {
            assert!(is_prime(&Integer::from(p as i64)), "{} should be prime", p);
        }
    }

    #[test]
    fn test_small_composites() {
        let composites = [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25];
        for &n in &composites {
            assert!(!is_prime(&Integer::from(n as i64)), "{} should be composite", n);
        }
    }

    #[test]
    fn test_edge_cases() {
        assert!(!is_prime(&Integer::from(0i64)));
        assert!(!is_prime(&Integer::from(1i64)));
        assert!(is_prime(&Integer::from(2i64)));
    }

    #[test]
    fn test_large_primes() {
        // Mersenne prime M61 = 2^61 - 1
        let m61 = Integer::from(2305843009213693951i64);
        assert!(is_prime(&m61));

        // A 64-bit prime
        let large = Integer::from(18446744073709551557u64);
        assert!(is_prime(&large));
    }

    #[test]
    fn test_carmichael_numbers() {
        // Carmichael numbers fool Fermat's test but not Miller-Rabin
        let carmichael = [561, 1105, 1729, 2465, 2821, 6601, 8911];
        for &n in &carmichael {
            assert!(!is_prime(&Integer::from(n as i64)),
                    "Carmichael number {} should be detected as composite", n);
        }
    }
}
