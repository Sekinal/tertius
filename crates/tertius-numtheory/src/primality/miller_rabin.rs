//! Miller-Rabin primality test.
//!
//! This is a probabilistic primality test that is very fast and reliable.
//! With k random witnesses, the probability of a false positive is at most 4^(-k).
//!
//! For deterministic testing of 64-bit numbers, we use specific witness sets
//! that are known to be sufficient.

use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use tertius_integers::Integer;

/// Performs k rounds of Miller-Rabin with random witnesses.
///
/// Returns true if n is probably prime, false if definitely composite.
/// False positive probability is at most 4^(-k).
pub fn miller_rabin(n: &Integer, k: u32) -> bool {
    // Handle small cases
    if n <= &Integer::from(1i64) {
        return false;
    }
    if n == &Integer::from(2i64) || n == &Integer::from(3i64) {
        return true;
    }
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return false;
    }

    // Write n-1 = 2^r * d where d is odd
    let n_minus_1 = n.clone() - Integer::one();
    let (r, d) = factor_out_twos(&n_minus_1);

    // Use a seeded RNG for reproducibility in tests
    let mut rng = ChaCha8Rng::seed_from_u64(0xDEADBEEF);

    for _ in 0..k {
        // Pick random witness a in [2, n-2]
        let a = random_in_range(&mut rng, &Integer::from(2i64), &(n.clone() - Integer::from(2i64)));
        if !miller_rabin_witness(n, &a, &d, r) {
            return false;
        }
    }

    true
}

/// Deterministic Miller-Rabin for 64-bit numbers.
///
/// Uses carefully chosen witnesses that are known to correctly classify
/// all 64-bit numbers.
pub fn miller_rabin_deterministic(n: &Integer) -> bool {
    // Handle small cases
    if n <= &Integer::from(1i64) {
        return false;
    }
    if n == &Integer::from(2i64) {
        return true;
    }
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return false;
    }

    // Write n-1 = 2^r * d where d is odd
    let n_minus_1 = n.clone() - Integer::one();
    let (r, d) = factor_out_twos(&n_minus_1);

    // Witnesses that work for all 64-bit numbers
    // Based on: https://miller-rabin.appspot.com/
    // These 7 witnesses are sufficient for n < 3,317,044,064,679,887,385,961,981
    let witnesses: &[u64] = if n < &Integer::from(2047u64) {
        &[2]
    } else if n < &Integer::from(1373653u64) {
        &[2, 3]
    } else if n < &Integer::from(9080191u64) {
        &[31, 73]
    } else if n < &Integer::from(25326001u64) {
        &[2, 3, 5]
    } else if n < &Integer::from(3215031751u64) {
        &[2, 3, 5, 7]
    } else if n < &Integer::from(4759123141u64) {
        &[2, 7, 61]
    } else if n < &Integer::from(1122004669633u64) {
        &[2, 13, 23, 1662803]
    } else if n < &Integer::from(2152302898747u64) {
        &[2, 3, 5, 7, 11]
    } else if n < &Integer::from(3474749660383u64) {
        &[2, 3, 5, 7, 11, 13]
    } else if n < &Integer::from(341550071728321u64) {
        &[2, 3, 5, 7, 11, 13, 17]
    } else if n < &Integer::from(3825123056546413051u64) {
        &[2, 3, 5, 7, 11, 13, 17, 19, 23]
    } else {
        // For n >= 3825123056546413051, use these witnesses
        // which work for all 64-bit numbers
        &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    };

    for &a in witnesses {
        let a_int = Integer::from(a);
        if &a_int >= n {
            continue;
        }
        if !miller_rabin_witness(n, &a_int, &d, r) {
            return false;
        }
    }

    true
}

/// Tests a single Miller-Rabin witness.
///
/// Returns true if n passes the test with witness a, false if n is composite.
fn miller_rabin_witness(n: &Integer, a: &Integer, d: &Integer, r: u32) -> bool {
    // Compute a^d mod n
    let mut x = mod_pow(a, d, n);

    let n_minus_1 = n.clone() - Integer::one();

    if x.is_one() || x == n_minus_1 {
        return true;
    }

    // Square r-1 times
    for _ in 0..r - 1 {
        x = mod_pow(&x, &Integer::from(2i64), n);
        if x == n_minus_1 {
            return true;
        }
        if x.is_one() {
            return false;
        }
    }

    false
}

/// Factors out powers of 2: returns (r, d) such that n = 2^r * d and d is odd.
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

/// Computes base^exp mod m using binary exponentiation.
pub(crate) fn mod_pow(base: &Integer, exp: &Integer, m: &Integer) -> Integer {
    if m.is_one() {
        return Integer::zero();
    }

    let mut result = Integer::one();
    let mut base = base.clone() % m.clone();
    let mut exp = exp.clone();
    let two = Integer::from(2i64);

    while !exp.is_zero() {
        // If exp is odd, multiply result by base
        if !(exp.clone() % two.clone()).is_zero() {
            result = (result * base.clone()) % m.clone();
        }
        exp = exp / two.clone();
        base = (base.clone() * base.clone()) % m.clone();
    }

    result
}

/// Generates a random Integer in the range [low, high].
fn random_in_range(rng: &mut ChaCha8Rng, low: &Integer, high: &Integer) -> Integer {
    if low >= high {
        return low.clone();
    }

    let range = high.clone() - low.clone() + Integer::one();

    // For simplicity, generate random bytes and reduce mod range
    // This has slight bias for very large ranges, but is fine for primality testing
    let bit_len = range.bit_len();
    let byte_len = (bit_len + 7) / 8;

    loop {
        let mut bytes = vec![0u8; byte_len];
        rng.fill(&mut bytes[..]);

        // Convert to Integer
        let mut r = Integer::zero();
        for byte in bytes {
            r = r * Integer::from(256i64) + Integer::from(byte as i64);
        }

        // Reduce mod range
        r = r % range.clone();

        // Return low + r
        return low.clone() + r;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_pow() {
        // 2^10 mod 1000 = 1024 mod 1000 = 24
        assert_eq!(
            mod_pow(&Integer::from(2i64), &Integer::from(10i64), &Integer::from(1000i64)),
            Integer::from(24i64)
        );

        // 3^7 mod 13 = 2187 mod 13 = 3
        assert_eq!(
            mod_pow(&Integer::from(3i64), &Integer::from(7i64), &Integer::from(13i64)),
            Integer::from(3i64)
        );
    }

    #[test]
    fn test_factor_out_twos() {
        // 24 = 2^3 * 3
        let (r, d) = factor_out_twos(&Integer::from(24i64));
        assert_eq!(r, 3);
        assert_eq!(d, Integer::from(3i64));

        // 15 = 2^0 * 15
        let (r, d) = factor_out_twos(&Integer::from(15i64));
        assert_eq!(r, 0);
        assert_eq!(d, Integer::from(15i64));
    }

    #[test]
    fn test_miller_rabin_deterministic() {
        // Small primes
        assert!(miller_rabin_deterministic(&Integer::from(2i64)));
        assert!(miller_rabin_deterministic(&Integer::from(3i64)));
        assert!(miller_rabin_deterministic(&Integer::from(5i64)));
        assert!(miller_rabin_deterministic(&Integer::from(997i64)));

        // Small composites
        assert!(!miller_rabin_deterministic(&Integer::from(4i64)));
        assert!(!miller_rabin_deterministic(&Integer::from(9i64)));
        assert!(!miller_rabin_deterministic(&Integer::from(1001i64)));

        // Carmichael numbers
        assert!(!miller_rabin_deterministic(&Integer::from(561i64)));
        assert!(!miller_rabin_deterministic(&Integer::from(1729i64)));
    }

    #[test]
    fn test_miller_rabin_large() {
        // 2^61 - 1 (Mersenne prime M61)
        let m61 = Integer::from(2305843009213693951i64);
        assert!(miller_rabin_deterministic(&m61));

        // 10^9 + 7 (common prime)
        assert!(miller_rabin_deterministic(&Integer::from(1000000007i64)));

        // 10^9 + 9 is also prime
        assert!(miller_rabin_deterministic(&Integer::from(1000000009i64)));

        // A composite: 10^9 + 8 = 2^3 * 125000001
        assert!(!miller_rabin_deterministic(&Integer::from(1000000008i64)));
    }

    #[test]
    fn test_miller_rabin_probabilistic() {
        assert!(miller_rabin(&Integer::from(1000000007i64), 25));
        assert!(!miller_rabin(&Integer::from(1000000008i64), 25));
    }
}
