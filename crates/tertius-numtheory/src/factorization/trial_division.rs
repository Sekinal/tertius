//! Trial division factorization.
//!
//! Simple but effective for extracting small prime factors.
//! This is typically used as a first step before more sophisticated methods.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::Factorization;

/// Small primes for trial division (primes up to 1000).
pub const SMALL_PRIMES: [u64; 168] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997,
];

/// Performs complete factorization using trial division.
///
/// This is slow for large numbers but guaranteed to find all factors.
/// Use this only for small numbers (< 10^12).
pub fn trial_division_factor(n: &Integer) -> Factorization {
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

    // Try small primes first
    for &p in &SMALL_PRIMES {
        let p_int = Integer::from(p);
        let mut exp = 0u32;
        while (n.clone() % p_int.clone()).is_zero() {
            n = n / p_int.clone();
            exp += 1;
        }
        if exp > 0 {
            factorization.add_factor(p_int, exp);
        }
        if n.is_one() {
            return factorization;
        }
    }

    // Continue with 6k±1 pattern for larger factors
    let mut k = 167u64; // 997 = 6*166 + 1, so start at 167
    let sqrt_n = integer_sqrt(&n);

    loop {
        let candidate = Integer::from(6 * k - 1);
        if &candidate > &sqrt_n {
            break;
        }
        let mut exp = 0u32;
        while (n.clone() % candidate.clone()).is_zero() {
            n = n / candidate.clone();
            exp += 1;
        }
        if exp > 0 {
            factorization.add_factor(candidate, exp);
        }

        let candidate = Integer::from(6 * k + 1);
        if &candidate > &sqrt_n {
            break;
        }
        let mut exp = 0u32;
        while (n.clone() % candidate.clone()).is_zero() {
            n = n / candidate.clone();
            exp += 1;
        }
        if exp > 0 {
            factorization.add_factor(candidate, exp);
        }

        k += 1;

        // Update sqrt_n after divisions
        if exp > 0 {
            let new_sqrt = integer_sqrt(&n);
            if new_sqrt < sqrt_n {
                // Continue checking
            }
        }
    }

    // If n > 1, then n is a prime factor
    if !n.is_one() {
        factorization.add_factor(n, 1);
    }

    factorization
}

/// Extracts small prime factors up to a bound.
///
/// Returns (remaining, factorization) where remaining is the unfactored part.
/// This is useful as a preprocessing step for more sophisticated methods.
pub fn extract_small_factors(n: &Integer, bound: u64) -> (Integer, Factorization) {
    if n.is_zero() {
        return (Integer::zero(), Factorization::zero());
    }

    let mut factorization = Factorization::one();
    let mut remaining = n.clone();

    if remaining.is_negative() {
        factorization.sign = -1;
        remaining = -remaining;
    }

    if remaining.is_one() {
        return (Integer::one(), factorization);
    }

    // Try small primes up to bound
    for &p in &SMALL_PRIMES {
        if p > bound {
            break;
        }
        let p_int = Integer::from(p);
        let mut exp = 0u32;
        while (remaining.clone() % p_int.clone()).is_zero() {
            remaining = remaining / p_int.clone();
            exp += 1;
        }
        if exp > 0 {
            factorization.add_factor(p_int, exp);
        }
        if remaining.is_one() {
            break;
        }
    }

    // Continue with 6k±1 if we haven't reached the bound
    if *SMALL_PRIMES.last().unwrap() < bound {
        let mut k = 167u64;
        loop {
            let candidate = 6 * k - 1;
            if candidate > bound {
                break;
            }
            let p_int = Integer::from(candidate);
            let mut exp = 0u32;
            while (remaining.clone() % p_int.clone()).is_zero() {
                remaining = remaining / p_int.clone();
                exp += 1;
            }
            if exp > 0 {
                factorization.add_factor(p_int, exp);
            }

            let candidate = 6 * k + 1;
            if candidate > bound {
                break;
            }
            let p_int = Integer::from(candidate);
            let mut exp = 0u32;
            while (remaining.clone() % p_int.clone()).is_zero() {
                remaining = remaining / p_int.clone();
                exp += 1;
            }
            if exp > 0 {
                factorization.add_factor(p_int, exp);
            }

            k += 1;
        }
    }

    (remaining, factorization)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_numbers() {
        assert!(trial_division_factor(&Integer::from(0i64)).value().is_zero());
        assert!(trial_division_factor(&Integer::from(1i64)).is_one());

        let f = trial_division_factor(&Integer::from(2i64));
        assert_eq!(f.to_vec(), vec![(Integer::from(2i64), 1)]);

        let f = trial_division_factor(&Integer::from(12i64));
        assert_eq!(
            f.to_vec(),
            vec![(Integer::from(2i64), 2), (Integer::from(3i64), 1)]
        );
    }

    #[test]
    fn test_prime() {
        let f = trial_division_factor(&Integer::from(997i64));
        assert_eq!(f.to_vec(), vec![(Integer::from(997i64), 1)]);
    }

    #[test]
    fn test_prime_power() {
        let f = trial_division_factor(&Integer::from(1024i64)); // 2^10
        assert_eq!(f.to_vec(), vec![(Integer::from(2i64), 10)]);
    }

    #[test]
    fn test_negative() {
        let f = trial_division_factor(&Integer::from(-30i64));
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
    fn test_large_composite() {
        // 1234567890 = 2 × 3^2 × 5 × 3607 × 3803
        let f = trial_division_factor(&Integer::from(1234567890i64));
        assert_eq!(f.value(), Integer::from(1234567890i64));
    }

    #[test]
    fn test_extract_small_factors() {
        // 2^10 × 3^5 × 10007 (10007 is prime > 1000)
        let n = Integer::from(2i64).pow(10) * Integer::from(3i64).pow(5) * Integer::from(10007i64);
        let (remaining, factors) = extract_small_factors(&n, 1000);
        assert_eq!(remaining, Integer::from(10007i64));
        assert_eq!(
            factors.to_vec(),
            vec![(Integer::from(2i64), 10), (Integer::from(3i64), 5)]
        );
    }
}
