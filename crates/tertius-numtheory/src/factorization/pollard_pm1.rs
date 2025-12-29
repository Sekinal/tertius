//! Pollard's p-1 factorization algorithm.
//!
//! This algorithm is effective when n has a prime factor p such that
//! p-1 is B-smooth (i.e., all prime factors of p-1 are ≤ B).
//!
//! The algorithm uses Fermat's little theorem: a^(p-1) ≡ 1 (mod p),
//! so a^M ≡ 1 (mod p) when (p-1) | M.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::trial_division::SMALL_PRIMES;

/// Pollard's p-1 algorithm with smoothness bound B.
///
/// Finds a factor of n if n has a prime factor p where p-1 is B-smooth.
/// The default bound B = 10^6 works well for most numbers.
pub fn pollard_pm1(n: &Integer) -> Option<Integer> {
    pollard_pm1_with_bound(n, 1_000_000)
}

/// Pollard's p-1 with custom smoothness bound.
pub fn pollard_pm1_with_bound(n: &Integer, b1: u64) -> Option<Integer> {
    if n.is_one() {
        return None;
    }

    // Check if n is even
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return Some(Integer::from(2i64));
    }

    // Start with a = 2 (could try other bases if this fails)
    let mut a = Integer::from(2i64);

    // Stage 1: Compute a^M mod n where M = lcm(1, 2, ..., B1)
    // We do this by computing a^(p^e) for each prime power p^e ≤ B1

    // First, handle small primes from our table
    for &p in &SMALL_PRIMES {
        if p > b1 {
            break;
        }

        // Find the largest e such that p^e ≤ B1
        let mut pe = p;
        while pe <= b1 / p {
            pe *= p;
        }

        // a = a^pe mod n
        a = mod_pow(&a, &Integer::from(pe), n);

        // Check if we found a factor
        let g = gcd(&(a.clone() - Integer::one()), n);
        if !g.is_one() && &g != n {
            return Some(g);
        }
    }

    // Continue with primes beyond our small prime table
    let start = *SMALL_PRIMES.last().unwrap() + 2;
    let mut p = start;
    while p <= b1 {
        if is_prime_simple(p) {
            // Find largest e such that p^e ≤ B1
            let mut pe = p;
            while pe <= b1 / p {
                pe *= p;
            }

            a = mod_pow(&a, &Integer::from(pe), n);

            // Periodically check for factors
            if p % 1000 == 1 {
                let g = gcd(&(a.clone() - Integer::one()), n);
                if !g.is_one() && &g != n {
                    return Some(g);
                }
            }
        }
        p += 2; // Skip even numbers
    }

    // Final GCD check
    let g = gcd(&(a.clone() - Integer::one()), n);
    if !g.is_one() && &g != n {
        return Some(g);
    }

    None
}

/// Two-stage Pollard p-1 with bounds B1 and B2.
///
/// Stage 1 handles prime powers up to B1.
/// Stage 2 handles primes in the range (B1, B2].
///
/// This is more efficient for finding factors with one large prime factor in p-1.
pub fn pollard_pm1_two_stage(n: &Integer, b1: u64, b2: u64) -> Option<Integer> {
    if n.is_one() {
        return None;
    }

    // Check if n is even
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return Some(Integer::from(2i64));
    }

    let mut a = Integer::from(2i64);

    // Stage 1: Same as single-stage
    for &p in &SMALL_PRIMES {
        if p > b1 {
            break;
        }
        let mut pe = p;
        while pe <= b1 / p {
            pe *= p;
        }
        a = mod_pow(&a, &Integer::from(pe), n);
    }

    // Check after stage 1
    let g = gcd(&(a.clone() - Integer::one()), n);
    if !g.is_one() && &g != n {
        return Some(g);
    }
    if &g == n {
        return None; // Stage 1 failed completely
    }

    // Stage 2: Handle primes in (B1, B2] one at a time
    // Use baby-step giant-step approach for efficiency
    let b = a.clone(); // b = a after stage 1

    // Precompute b^2, b^4, b^6, ... up to b^(gap)
    // where gap is the maximum prime gap in our range
    let max_gap = 200u64; // Conservative upper bound
    let mut powers = vec![Integer::one()]; // powers[i] = b^(2i)
    let mut b_sq = (b.clone() * b.clone()) % n.clone();
    for _ in 0..max_gap / 2 {
        powers.push(powers.last().unwrap().clone() * b_sq.clone() % n.clone());
    }

    // Process primes in (B1, B2]
    let mut current = b.clone();
    let mut last_prime = b1;
    let mut prod = Integer::one();

    // Simple sieve for primes in range
    let sieve_start = if b1 % 2 == 0 { b1 + 1 } else { b1 + 2 };
    for p in (sieve_start..=b2).step_by(2) {
        if !is_prime_simple(p) {
            continue;
        }

        let gap = (p - last_prime) as usize;
        if gap / 2 < powers.len() {
            // current = current * b^gap = current * (b^2)^(gap/2)
            current = (current * powers[gap / 2].clone()) % n.clone();
        } else {
            // Gap too large, compute directly
            current = mod_pow(&current, &Integer::from(p - last_prime), n);
        }

        prod = (prod * (current.clone() - Integer::one())) % n.clone();

        // Periodically check GCD
        if p % 10000 < 2 {
            let g = gcd(&prod, n);
            if !g.is_one() && &g != n {
                return Some(g);
            }
            prod = Integer::one();
        }

        last_prime = p;
    }

    // Final GCD
    let g = gcd(&prod, n);
    if !g.is_one() && &g != n {
        return Some(g);
    }

    None
}

/// Simple primality check for small numbers.
fn is_prime_simple(n: u64) -> bool {
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

/// Modular exponentiation: base^exp mod m.
fn mod_pow(base: &Integer, exp: &Integer, m: &Integer) -> Integer {
    if m.is_one() {
        return Integer::zero();
    }

    let mut result = Integer::one();
    let mut base = base.clone() % m.clone();
    let mut exp = exp.clone();
    let two = Integer::from(2i64);

    while !exp.is_zero() {
        if !(exp.clone() % two.clone()).is_zero() {
            result = (result * base.clone()) % m.clone();
        }
        exp = exp / two.clone();
        base = (base.clone() * base.clone()) % m.clone();
    }

    result
}

/// GCD using Euclidean algorithm.
fn gcd(a: &Integer, b: &Integer) -> Integer {
    a.gcd(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pollard_pm1_smooth() {
        // p = 1009, p-1 = 1008 = 2^4 × 3^2 × 7 (16-smooth)
        // q = 10007, q-1 = 10006 = 2 × 5003 (5003-smooth, 5003 is prime)
        // n = p × q should be factorable with B1 >= 5003

        // Simpler test: p = 1009 (1008 = 2^4 × 3^2 × 7)
        let n = Integer::from(1009i64) * Integer::from(1013i64);
        let factor = pollard_pm1_with_bound(&n, 100);
        // 1008 = 16 × 63 = 16 × 9 × 7, all factors ≤ 16
        // So B1 = 100 should work
        assert!(factor.is_some());
    }

    #[test]
    fn test_pollard_pm1_two_stage() {
        // Test with a number that needs two stages
        // p = 10007, p-1 = 2 × 5003 (5003 is prime)
        // With B1 = 100, B2 = 10000, should find factor in stage 2
        let n = Integer::from(10007i64) * Integer::from(10009i64);
        let factor = pollard_pm1_two_stage(&n, 100, 10000);
        assert!(factor.is_some());
    }
}
