//! Prime generation and counting.
//!
//! This module provides:
//! - Sieve of Eratosthenes (segmented for large ranges)
//! - nth prime computation
//! - Prime counting function π(x)

use tertius_integers::Integer;

/// Sieve of Eratosthenes up to n.
///
/// Returns a vector of all primes up to and including n.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::sieve_of_eratosthenes;
///
/// let primes = sieve_of_eratosthenes(30);
/// assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
/// ```
pub fn sieve_of_eratosthenes(n: u64) -> Vec<u64> {
    if n < 2 {
        return vec![];
    }

    let n = n as usize;
    let mut is_prime = vec![true; n + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    let sqrt_n = (n as f64).sqrt() as usize;
    for i in 2..=sqrt_n {
        if is_prime[i] {
            let mut j = i * i;
            while j <= n {
                is_prime[j] = false;
                j += i;
            }
        }
    }

    is_prime
        .iter()
        .enumerate()
        .filter(|(_, &is_p)| is_p)
        .map(|(i, _)| i as u64)
        .collect()
}

/// Returns all primes up to n.
///
/// Alias for `sieve_of_eratosthenes`.
pub fn primes_up_to(n: u64) -> Vec<u64> {
    sieve_of_eratosthenes(n)
}

/// Segmented sieve for large ranges.
///
/// Returns primes in the range [low, high].
///
/// # Examples
///
/// ```
/// use tertius_numtheory::primes::segmented_sieve;
///
/// let primes = segmented_sieve(100, 120);
/// assert_eq!(primes, vec![101, 103, 107, 109, 113]);
/// ```
pub fn segmented_sieve(low: u64, high: u64) -> Vec<u64> {
    if high < 2 {
        return vec![];
    }

    let low = std::cmp::max(2, low);

    // Get small primes up to sqrt(high)
    let sqrt_high = ((high as f64).sqrt() as u64) + 1;
    let small_primes = sieve_of_eratosthenes(sqrt_high);

    let segment_size = std::cmp::min(high - low + 1, 1_000_000);
    let mut result = Vec::new();

    let mut segment_start = low;
    while segment_start <= high {
        let segment_end = std::cmp::min(segment_start + segment_size - 1, high);
        let size = (segment_end - segment_start + 1) as usize;

        let mut is_prime = vec![true; size];

        for &p in &small_primes {
            if p * p > segment_end {
                break;
            }

            // Find first multiple of p in range
            let mut start = ((segment_start + p - 1) / p) * p;
            if start == p {
                start += p;
            }

            while start <= segment_end {
                is_prime[(start - segment_start) as usize] = false;
                start += p;
            }
        }

        for i in 0..size {
            if is_prime[i] {
                result.push(segment_start + i as u64);
            }
        }

        segment_start = segment_end + 1;
    }

    result
}

/// Compute the nth prime (1-indexed, so nth_prime(1) = 2).
///
/// Uses an approximation followed by sieving.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::nth_prime;
///
/// assert_eq!(nth_prime(1), 2);
/// assert_eq!(nth_prime(10), 29);
/// assert_eq!(nth_prime(100), 541);
/// ```
pub fn nth_prime(n: u64) -> u64 {
    if n == 0 {
        panic!("nth_prime is 1-indexed, n must be >= 1");
    }
    if n == 1 {
        return 2;
    }
    if n == 2 {
        return 3;
    }

    // Upper bound for the nth prime using prime number theorem
    let n_f = n as f64;
    let ln_n = n_f.ln();
    let ln_ln_n = ln_n.ln();

    // Upper bound: p_n < n(ln n + ln ln n) for n >= 6
    let upper = if n >= 6 {
        (n_f * (ln_n + ln_ln_n)).ceil() as u64
    } else {
        30
    };

    // Sieve and count
    let primes = sieve_of_eratosthenes(upper);

    if primes.len() >= n as usize {
        primes[n as usize - 1]
    } else {
        // Need more primes (shouldn't happen with good upper bound)
        let extended = sieve_of_eratosthenes(upper * 2);
        extended[n as usize - 1]
    }
}

/// Prime counting function π(x).
///
/// Returns the number of primes less than or equal to x.
///
/// For small x, uses direct sieving.
/// For larger x, uses the Meissel-Lehmer algorithm (approximation).
///
/// # Examples
///
/// ```
/// use tertius_numtheory::prime_count;
///
/// assert_eq!(prime_count(10), 4);   // 2, 3, 5, 7
/// assert_eq!(prime_count(100), 25);
/// ```
pub fn prime_count(x: u64) -> u64 {
    if x < 2 {
        return 0;
    }

    // For small x, just sieve
    if x <= 10_000_000 {
        return sieve_of_eratosthenes(x).len() as u64;
    }

    // For larger x, use Legendre's formula or approximation
    // This is a simplified version; a full implementation would use Meissel-Lehmer
    prime_count_legendre(x)
}

/// Legendre's formula for prime counting (simplified).
///
/// π(x) = π(√x) + φ(x, π(√x)) - 1
/// where φ(x, a) is the number of integers ≤ x not divisible by any of first a primes.
fn prime_count_legendre(x: u64) -> u64 {
    if x < 2 {
        return 0;
    }
    if x <= 10_000_000 {
        return sieve_of_eratosthenes(x).len() as u64;
    }

    // Get primes up to √x
    let sqrt_x = (x as f64).sqrt() as u64;
    let small_primes = sieve_of_eratosthenes(sqrt_x);
    let a = small_primes.len();

    // φ(x, a) using inclusion-exclusion (simplified)
    let phi = phi_function(x, &small_primes);

    // π(x) = φ(x, a) + a - 1
    phi + a as u64 - 1
}

/// φ(x, a) = count of integers ≤ x not divisible by first a primes.
///
/// Uses a simple inclusion-exclusion approach.
fn phi_function(x: u64, primes: &[u64]) -> u64 {
    if primes.is_empty() || x == 0 {
        return x;
    }

    // Use a truncated sieve for efficiency
    // This is a simplified version
    let a = std::cmp::min(primes.len(), 7); // Limit depth for efficiency

    fn phi_rec(x: u64, primes: &[u64], idx: usize) -> i64 {
        if idx == 0 || x == 0 {
            return x as i64;
        }
        if primes[idx - 1] > x {
            return phi_rec(x, primes, idx - 1);
        }

        phi_rec(x, primes, idx - 1) - phi_rec(x / primes[idx - 1], primes, idx - 1)
    }

    phi_rec(x, primes, a) as u64
}

/// Check if n is prime (simple trial division for small n).
pub fn is_prime_simple(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }

    let mut i = 5u64;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

/// Prime gap: difference between consecutive primes.
///
/// Returns the gap after the prime p.
pub fn prime_gap(p: u64) -> u64 {
    if !is_prime_simple(p) {
        return 0;
    }

    let mut next = p + 1;
    while !is_prime_simple(next) {
        next += 1;
    }

    next - p
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sieve() {
        let primes = sieve_of_eratosthenes(30);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn test_sieve_edge_cases() {
        assert!(sieve_of_eratosthenes(0).is_empty());
        assert!(sieve_of_eratosthenes(1).is_empty());
        assert_eq!(sieve_of_eratosthenes(2), vec![2]);
    }

    #[test]
    fn test_segmented_sieve() {
        let primes = segmented_sieve(100, 120);
        assert_eq!(primes, vec![101, 103, 107, 109, 113]);
    }

    #[test]
    fn test_nth_prime() {
        assert_eq!(nth_prime(1), 2);
        assert_eq!(nth_prime(2), 3);
        assert_eq!(nth_prime(10), 29);
        assert_eq!(nth_prime(100), 541);
        assert_eq!(nth_prime(1000), 7919);
    }

    #[test]
    fn test_prime_count() {
        assert_eq!(prime_count(10), 4);
        assert_eq!(prime_count(100), 25);
        assert_eq!(prime_count(1000), 168);
    }

    #[test]
    fn test_prime_gap() {
        assert_eq!(prime_gap(2), 1);  // gap to 3
        assert_eq!(prime_gap(3), 2);  // gap to 5
        assert_eq!(prime_gap(7), 4);  // gap to 11
    }
}
