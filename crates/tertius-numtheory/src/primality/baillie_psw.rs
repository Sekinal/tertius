//! Baillie-PSW primality test.
//!
//! The Baillie-PSW test combines Miller-Rabin with base 2 and a strong Lucas test.
//! There are no known pseudoprimes for this test, making it effectively deterministic
//! for all practical purposes.
//!
//! For 64-bit numbers, this test is provably correct (no counterexamples exist).
//! For larger numbers, it is probabilistic but has never been wrong.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::miller_rabin::{miller_rabin_deterministic, mod_pow};
use super::lucas::strong_lucas_test;
use super::trial_division::quick_check;
use super::PrimalityResult;

/// Baillie-PSW primality test.
///
/// This is the recommended primality test for general use.
/// It is:
/// - Deterministic for n < 2^64 (no known pseudoprimes)
/// - Extremely reliable for larger numbers (no known pseudoprimes at all)
/// - Fast: O(log³n) operations
///
/// # Examples
///
/// ```
/// use tertius_numtheory::primality::baillie_psw;
/// use tertius_integers::Integer;
///
/// // Test a prime
/// let result = baillie_psw(&Integer::from(1000000007i64));
/// assert!(result.is_prime());
///
/// // Test a composite
/// let result = baillie_psw(&Integer::from(561i64)); // Carmichael number
/// assert!(!result.is_prime());
/// ```
pub fn baillie_psw(n: &Integer) -> PrimalityResult {
    // Quick check for small factors and special cases
    if let Some(result) = quick_check(n) {
        return result;
    }

    // Check if n is a perfect square (perfect squares > 1 are composite)
    if is_perfect_square(n) {
        return PrimalityResult::Composite;
    }

    // Step 1: Miller-Rabin with base 2
    // This is a quick filter that catches most composites
    if !miller_rabin_base_2(n) {
        return PrimalityResult::Composite;
    }

    // Step 2: Strong Lucas test
    // This catches composites that passed Miller-Rabin base 2
    if !strong_lucas_test(n) {
        return PrimalityResult::Composite;
    }

    // If both tests pass, n is almost certainly prime
    // For n < 2^64, this is provably correct
    // For larger n, there are no known counterexamples
    if n.bit_len() <= 64 {
        PrimalityResult::Prime
    } else {
        PrimalityResult::ProbablyPrime
    }
}

/// Miller-Rabin test with base 2 only.
fn miller_rabin_base_2(n: &Integer) -> bool {
    if n <= &Integer::from(2i64) {
        return n == &Integer::from(2i64);
    }
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return false;
    }

    // Write n-1 = 2^r * d where d is odd
    let n_minus_1 = n.clone() - Integer::one();
    let (r, d) = factor_out_twos(&n_minus_1);

    // Test with base 2
    let mut x = mod_pow(&Integer::from(2i64), &d, n);

    if x.is_one() || x == n_minus_1 {
        return true;
    }

    for _ in 0..r - 1 {
        x = (x.clone() * x.clone()) % n.clone();
        if x == n_minus_1 {
            return true;
        }
        if x.is_one() {
            return false;
        }
    }

    false
}

/// Check if n is a perfect square.
fn is_perfect_square(n: &Integer) -> bool {
    if n.is_zero() || n.is_one() {
        return true;
    }
    if n.is_negative() {
        return false;
    }

    let sqrt = integer_sqrt(n);
    &(sqrt.clone() * sqrt) == n
}

/// Integer square root using Newton's method.
fn integer_sqrt(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }

    let mut x = n.clone();
    let mut y = (n.clone() + Integer::one()) / Integer::from(2i64);

    while y < x {
        x = y.clone();
        y = (x.clone() + n.clone() / x.clone()) / Integer::from(2i64);
    }

    x
}

/// Factor out powers of 2.
fn factor_out_twos(n: &Integer) -> (u32, Integer) {
    let mut d = n.clone();
    let mut r = 0u32;
    let two = Integer::from(2i64);

    while (d.clone() % two.clone()).is_zero() {
        d = d / two.clone();
        r += 1;
    }

    (r, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_primes() {
        let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
        for &p in &primes {
            assert!(
                baillie_psw(&Integer::from(p as i64)).is_prime(),
                "{} should be prime",
                p
            );
        }
    }

    #[test]
    fn test_composites() {
        let composites = [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26];
        for &n in &composites {
            assert!(
                !baillie_psw(&Integer::from(n as i64)).is_prime(),
                "{} should be composite",
                n
            );
        }
    }

    #[test]
    fn test_carmichael_numbers() {
        // Carmichael numbers are pseudoprimes for Fermat's test
        // but Baillie-PSW correctly identifies them as composite
        let carmichael = [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841];
        for &n in &carmichael {
            assert!(
                !baillie_psw(&Integer::from(n as i64)).is_prime(),
                "Carmichael number {} should be composite",
                n
            );
        }
    }

    #[test]
    fn test_perfect_squares() {
        // Perfect squares > 1 are composite
        assert!(!baillie_psw(&Integer::from(4i64)).is_prime());
        assert!(!baillie_psw(&Integer::from(9i64)).is_prime());
        assert!(!baillie_psw(&Integer::from(16i64)).is_prime());
        assert!(!baillie_psw(&Integer::from(10000i64)).is_prime());
    }

    #[test]
    fn test_large_primes() {
        // 2^61 - 1 (Mersenne prime)
        let m61 = Integer::from(2305843009213693951i64);
        assert!(baillie_psw(&m61).is_prime());

        // 10^9 + 7
        assert!(baillie_psw(&Integer::from(1000000007i64)).is_prime());

        // A large 64-bit prime
        let large = Integer::from(18446744073709551557u64);
        assert!(baillie_psw(&large).is_prime());
    }

    #[test]
    fn test_lucas_pseudoprimes() {
        // These are pseudoprimes for some Lucas tests but not for BPSW
        // 5459 is a Fibonacci pseudoprime
        assert!(!baillie_psw(&Integer::from(5459i64)).is_prime());
    }

    #[test]
    fn test_is_perfect_square() {
        assert!(is_perfect_square(&Integer::from(0i64)));
        assert!(is_perfect_square(&Integer::from(1i64)));
        assert!(is_perfect_square(&Integer::from(4i64)));
        assert!(is_perfect_square(&Integer::from(9i64)));
        assert!(is_perfect_square(&Integer::from(100i64)));

        assert!(!is_perfect_square(&Integer::from(2i64)));
        assert!(!is_perfect_square(&Integer::from(3i64)));
        assert!(!is_perfect_square(&Integer::from(5i64)));
        assert!(!is_perfect_square(&Integer::from(99i64)));
    }
}
