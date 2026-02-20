//! Factorization dispatcher.
//!
//! Automatically selects the best algorithm based on input size.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::ecm::ecm_factor;
use super::pollard_pm1::pollard_pm1;
use super::pollard_rho::pollard_rho_brent;
use super::siqs::siqs_factor;
use super::trial_division::{extract_small_factors, trial_division_factor};
use super::Factorization;
use crate::primality::is_prime;

/// Factor an integer using the best available algorithm.
///
/// The algorithm is automatically selected based on the size of n:
///
/// | Size | Algorithm |
/// |------|-----------|
/// | < 10^6 | Trial division |
/// | 10^6 - 10^15 | Pollard's rho |
/// | 10^15 - 10^25 | ECM |
/// | 10^25 - 10^100 | SIQS |
///
/// # Examples
///
/// ```
/// use tertius_numtheory::factor;
/// use tertius_integers::Integer;
///
/// let n = Integer::from(1234567890i64);
/// let factors = factor(&n);
/// assert_eq!(factors.value(), n);
/// ```
pub fn factor(n: &Integer) -> Factorization {
    // Handle special cases
    if n.is_zero() {
        return Factorization::zero();
    }

    let mut factorization = Factorization::one();

    // Handle sign
    let mut n = n.clone();
    if n.is_negative() {
        factorization.sign = -1;
        n = -n;
    }

    if n.is_one() {
        return factorization;
    }

    // Factor the number
    factor_recursive(&n, &mut factorization);

    factorization
}

/// Recursive factorization helper.
fn factor_recursive(n: &Integer, factorization: &mut Factorization) {
    if n.is_one() {
        return;
    }

    // Check if n is prime
    if is_prime(n) {
        factorization.add_factor(n.clone(), 1);
        return;
    }

    // Extract small factors first
    let (remaining, small_factors) = extract_small_factors(n, 10000);
    for (p, e) in small_factors.iter() {
        factorization.add_factor(p.clone(), *e);
    }

    if remaining.is_one() {
        return;
    }

    if is_prime(&remaining) {
        factorization.add_factor(remaining, 1);
        return;
    }

    // Find a non-trivial factor
    if let Some(factor) = find_factor(&remaining) {
        // Recursively factor both parts
        factor_recursive(&factor, factorization);

        let cofactor = remaining / factor;
        factor_recursive(&cofactor, factorization);
    } else {
        // Fallback: treat as prime (shouldn't happen)
        factorization.add_factor(remaining, 1);
    }
}

/// Find a non-trivial factor using the best available algorithm.
fn find_factor(n: &Integer) -> Option<Integer> {
    let bits = n.bit_len();

    // Small numbers: use Pollard's rho
    if bits <= 50 {
        return pollard_rho_brent(n);
    }

    // Medium numbers: try Pollard's rho first, then p-1
    if bits <= 80 {
        if let Some(f) = pollard_rho_brent(n) {
            return Some(f);
        }
        if let Some(f) = pollard_pm1(n) {
            return Some(f);
        }
        return ecm_factor(n);
    }

    // Large numbers: ECM for smooth factors, then SIQS
    if bits <= 200 {
        // Try quick Pollard methods first
        if let Some(f) = pollard_rho_brent(n) {
            return Some(f);
        }
        if let Some(f) = pollard_pm1(n) {
            return Some(f);
        }

        // ECM for medium-sized factors
        if let Some(f) = ecm_factor(n) {
            return Some(f);
        }

        // Fall back to SIQS
        return siqs_factor(n);
    }

    // Very large numbers: SIQS
    if let Some(f) = pollard_rho_brent(n) {
        return Some(f);
    }
    if let Some(f) = ecm_factor(n) {
        return Some(f);
    }
    siqs_factor(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factor_small() {
        let f = factor(&Integer::from(12i64));
        assert_eq!(f.value(), Integer::from(12i64));
        assert_eq!(
            f.to_vec(),
            vec![(Integer::from(2i64), 2), (Integer::from(3i64), 1)]
        );
    }

    #[test]
    fn test_factor_prime() {
        let f = factor(&Integer::from(1000000007i64));
        assert_eq!(f.to_vec(), vec![(Integer::from(1000000007i64), 1)]);
    }

    #[test]
    fn test_factor_semiprime() {
        let n = Integer::from(1000003i64) * Integer::from(1000033i64);
        let f = factor(&n);
        assert_eq!(f.value(), n);
        assert_eq!(f.num_distinct_factors(), 2);
    }

    #[test]
    fn test_factor_power_of_two() {
        let f = factor(&Integer::from(1024i64));
        assert_eq!(f.to_vec(), vec![(Integer::from(2i64), 10)]);
    }

    #[test]
    fn test_factor_negative() {
        let f = factor(&Integer::from(-30i64));
        assert_eq!(f.sign, -1);
        assert_eq!(
            f.to_vec(),
            vec![
                (Integer::from(2i64), 1),
                (Integer::from(3i64), 1),
                (Integer::from(5i64), 1)
            ]
        );
    }

    #[test]
    fn test_factor_large() {
        // A product of two 10-digit primes
        let n = Integer::from(1000000007i64) * Integer::from(1000000009i64);
        let f = factor(&n);
        assert_eq!(f.value(), n);
    }
}
