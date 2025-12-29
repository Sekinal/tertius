//! Trial division primality testing.
//!
//! Simple but effective for small numbers. Tests divisibility by 2, 3,
//! and numbers of the form 6k±1 up to √n.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::PrimalityResult;

/// Small primes for quick divisibility checks.
const SMALL_PRIMES: [u64; 54] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251,
];

/// Product of first few primes for GCD-based divisibility check.
const PRIMORIAL_SMALL: u64 = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23;

/// Quick check for small factors and special cases.
///
/// Returns `Some(result)` if the primality can be determined quickly,
/// `None` if more sophisticated tests are needed.
pub fn quick_check(n: &Integer) -> Option<PrimalityResult> {
    // Handle small cases
    if n.is_zero() || n.is_one() {
        return Some(PrimalityResult::Composite);
    }

    // Handle negative numbers
    if n.is_negative() {
        return Some(PrimalityResult::Composite);
    }

    // Check if it's a small prime
    if let Some(small) = n.to_i64() {
        if small <= 1 {
            return Some(PrimalityResult::Composite);
        }
        if small == 2 || small == 3 {
            return Some(PrimalityResult::Prime);
        }
        if small < 251 {
            for &p in &SMALL_PRIMES {
                if small == p as i64 {
                    return Some(PrimalityResult::Prime);
                }
                if small < p as i64 {
                    return Some(PrimalityResult::Composite);
                }
            }
        }
    }

    // Check divisibility by small primes
    for &p in &SMALL_PRIMES {
        let p_int = Integer::from(p);
        if n == &p_int {
            return Some(PrimalityResult::Prime);
        }
        let rem = n.clone() % p_int;
        if rem.is_zero() {
            return Some(PrimalityResult::Composite);
        }
    }

    None
}

/// Full trial division up to √n or the specified limit.
///
/// This is exact but slow for large numbers.
/// Use this only when you need certainty and n is relatively small.
pub fn trial_division(n: &Integer) -> PrimalityResult {
    // Handle quick cases
    if let Some(result) = quick_check(n) {
        return result;
    }

    // For numbers that passed the quick check, we need to continue
    // trial division beyond the small primes.
    let sqrt_n = integer_sqrt(n);

    // Continue with 6k±1 pattern starting after our small prime list
    let mut k = 42u64; // 6*42 = 252, just past 251
    loop {
        let candidate1 = Integer::from(6 * k - 1);
        if &candidate1 > &sqrt_n {
            break;
        }
        if (n.clone() % candidate1).is_zero() {
            return PrimalityResult::Composite;
        }

        let candidate2 = Integer::from(6 * k + 1);
        if &candidate2 > &sqrt_n {
            break;
        }
        if (n.clone() % candidate2).is_zero() {
            return PrimalityResult::Composite;
        }

        k += 1;
    }

    PrimalityResult::Prime
}

/// Computes the integer square root of n (floor(√n)).
fn integer_sqrt(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }

    // Newton's method
    let mut x = n.clone();
    let mut y = (n.clone() + Integer::one()) / Integer::from(2i64);

    while y < x {
        x = y.clone();
        y = (x.clone() + n.clone() / x.clone()) / Integer::from(2i64);
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_check_small() {
        assert_eq!(quick_check(&Integer::from(0i64)), Some(PrimalityResult::Composite));
        assert_eq!(quick_check(&Integer::from(1i64)), Some(PrimalityResult::Composite));
        assert_eq!(quick_check(&Integer::from(2i64)), Some(PrimalityResult::Prime));
        assert_eq!(quick_check(&Integer::from(3i64)), Some(PrimalityResult::Prime));
        assert_eq!(quick_check(&Integer::from(4i64)), Some(PrimalityResult::Composite));
        assert_eq!(quick_check(&Integer::from(5i64)), Some(PrimalityResult::Prime));
    }

    #[test]
    fn test_quick_check_divisibility() {
        // 1001 = 7 * 11 * 13, should be caught by quick check
        assert_eq!(quick_check(&Integer::from(1001i64)), Some(PrimalityResult::Composite));
    }

    #[test]
    fn test_trial_division() {
        assert_eq!(trial_division(&Integer::from(997i64)), PrimalityResult::Prime);
        assert_eq!(trial_division(&Integer::from(1009i64)), PrimalityResult::Prime);
        assert_eq!(trial_division(&Integer::from(1000i64)), PrimalityResult::Composite);
    }

    #[test]
    fn test_integer_sqrt() {
        assert_eq!(integer_sqrt(&Integer::from(0i64)), Integer::from(0i64));
        assert_eq!(integer_sqrt(&Integer::from(1i64)), Integer::from(1i64));
        assert_eq!(integer_sqrt(&Integer::from(4i64)), Integer::from(2i64));
        assert_eq!(integer_sqrt(&Integer::from(9i64)), Integer::from(3i64));
        assert_eq!(integer_sqrt(&Integer::from(10i64)), Integer::from(3i64));
        assert_eq!(integer_sqrt(&Integer::from(100i64)), Integer::from(10i64));
    }
}
