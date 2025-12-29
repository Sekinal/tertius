//! Integer factorization algorithms.
//!
//! This module provides multiple factorization algorithms that are automatically
//! selected based on the size of the input:
//!
//! | Size | Algorithm | Performance |
//! |------|-----------|-------------|
//! | < 10^6 | Trial division | Instant |
//! | 10^6 - 10^15 | Pollard's rho | < 1 second |
//! | 10^15 - 10^25 | ECM | Seconds to minutes |
//! | 10^25 - 10^100 | SIQS | Minutes to hours |
//!
//! The main entry point is [`factor`] which automatically selects the best algorithm.

mod trial_division;
mod pollard_rho;
mod pollard_pm1;
mod ecm;
mod siqs;
mod dispatcher;

pub use trial_division::trial_division_factor;
pub use pollard_rho::{pollard_rho, pollard_rho_brent};
pub use pollard_pm1::pollard_pm1;
pub use ecm::ecm_factor;
pub use siqs::siqs_factor;
pub use dispatcher::factor;

use std::collections::BTreeMap;
use tertius_integers::Integer;

/// A complete factorization of an integer.
///
/// Represents n = ±∏ pᵢ^eᵢ where pᵢ are distinct primes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Factorization {
    /// The sign of the original number (-1 or 1).
    pub sign: i8,
    /// Map from prime factors to their exponents.
    pub factors: BTreeMap<Integer, u32>,
}

impl Factorization {
    /// Creates a new empty factorization (representing 1).
    pub fn one() -> Self {
        Self {
            sign: 1,
            factors: BTreeMap::new(),
        }
    }

    /// Creates a factorization representing zero.
    pub fn zero() -> Self {
        Self {
            sign: 0,
            factors: BTreeMap::new(),
        }
    }

    /// Creates a factorization from a prime.
    pub fn from_prime(p: Integer) -> Self {
        let mut factors = BTreeMap::new();
        factors.insert(p, 1);
        Self { sign: 1, factors }
    }

    /// Adds a prime factor with multiplicity.
    pub fn add_factor(&mut self, p: Integer, exp: u32) {
        if exp > 0 {
            *self.factors.entry(p).or_insert(0) += exp;
        }
    }

    /// Multiplies this factorization by p^exp.
    pub fn mul_prime(&mut self, p: Integer, exp: u32) {
        self.add_factor(p, exp);
    }

    /// Returns the number represented by this factorization.
    pub fn value(&self) -> Integer {
        if self.sign == 0 {
            return Integer::from(0i64);
        }

        let mut result = Integer::from(1i64);
        for (p, &exp) in &self.factors {
            result = result * p.pow(exp);
        }

        if self.sign < 0 {
            -result
        } else {
            result
        }
    }

    /// Returns the number of distinct prime factors (ω(n)).
    pub fn num_distinct_factors(&self) -> usize {
        self.factors.len()
    }

    /// Returns the total number of prime factors with multiplicity (Ω(n)).
    pub fn num_factors_with_multiplicity(&self) -> u32 {
        self.factors.values().sum()
    }

    /// Returns an iterator over (prime, exponent) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&Integer, &u32)> {
        self.factors.iter()
    }

    /// Converts to a vector of (prime, exponent) pairs.
    pub fn to_vec(&self) -> Vec<(Integer, u32)> {
        self.factors.iter().map(|(p, &e)| (p.clone(), e)).collect()
    }

    /// Checks if the factorization represents 1.
    pub fn is_one(&self) -> bool {
        self.sign == 1 && self.factors.is_empty()
    }

    /// Checks if the number is a prime power.
    pub fn is_prime_power(&self) -> bool {
        self.factors.len() == 1
    }

    /// Checks if the number is squarefree.
    pub fn is_squarefree(&self) -> bool {
        self.factors.values().all(|&e| e == 1)
    }
}

impl std::fmt::Display for Factorization {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.sign == 0 {
            return write!(f, "0");
        }
        if self.factors.is_empty() {
            return write!(f, "{}", if self.sign < 0 { "-1" } else { "1" });
        }

        if self.sign < 0 {
            write!(f, "-1 × ")?;
        }

        let mut first = true;
        for (p, &exp) in &self.factors {
            if !first {
                write!(f, " × ")?;
            }
            first = false;
            if exp == 1 {
                write!(f, "{}", p)?;
            } else {
                write!(f, "{}^{}", p, exp)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorization_display() {
        let mut f = Factorization::one();
        f.add_factor(Integer::from(2i64), 3);
        f.add_factor(Integer::from(3i64), 2);
        f.add_factor(Integer::from(5i64), 1);
        // 2^3 × 3^2 × 5 = 360
        assert_eq!(f.value(), Integer::from(360i64));
        assert_eq!(f.to_string(), "2^3 × 3^2 × 5");
    }

    #[test]
    fn test_factorization_properties() {
        let mut f = Factorization::one();
        f.add_factor(Integer::from(2i64), 1);
        f.add_factor(Integer::from(3i64), 1);
        f.add_factor(Integer::from(5i64), 1);
        // 30 = 2 × 3 × 5
        assert_eq!(f.num_distinct_factors(), 3);
        assert_eq!(f.num_factors_with_multiplicity(), 3);
        assert!(f.is_squarefree());
        assert!(!f.is_prime_power());
    }

    #[test]
    fn test_prime_power() {
        let mut f = Factorization::one();
        f.add_factor(Integer::from(7i64), 4);
        // 7^4 = 2401
        assert!(f.is_prime_power());
        assert!(!f.is_squarefree());
        assert_eq!(f.value(), Integer::from(2401i64));
    }
}
