//! Pollard's rho factorization algorithm.
//!
//! This is a probabilistic factorization algorithm with expected running time
//! O(n^1/4) for finding a non-trivial factor of n.
//!
//! We implement both the classic Floyd's cycle detection and Brent's
//! improvement, which is about 24% faster on average.

use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use tertius_integers::Integer;

/// Pollard's rho algorithm with Floyd's cycle detection.
///
/// Attempts to find a non-trivial factor of n.
/// Returns `Some(factor)` if successful, `None` if the algorithm fails
/// (which can happen for prime numbers or with unlucky parameters).
pub fn pollard_rho(n: &Integer) -> Option<Integer> {
    pollard_rho_with_params(n, Integer::from(2i64), Integer::from(1i64))
}

/// Pollard's rho with customizable starting point and constant.
///
/// Uses the iteration f(x) = x² + c (mod n).
pub fn pollard_rho_with_params(n: &Integer, x0: Integer, c: Integer) -> Option<Integer> {
    if n.is_one() {
        return None;
    }

    // Check if n is even
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return Some(Integer::from(2i64));
    }

    let mut x = x0.clone();
    let mut y = x0;
    let mut d = Integer::one();

    let f = |x: Integer| -> Integer { ((x.clone() * x.clone()) + c.clone()) % n.clone() };

    // Floyd's cycle detection: tortoise and hare
    while d.is_one() {
        x = f(x); // Tortoise moves one step
        y = f(f(y)); // Hare moves two steps

        // d = gcd(|x - y|, n)
        let diff = if x >= y {
            x.clone() - y.clone()
        } else {
            y.clone() - x.clone()
        };
        d = gcd(&diff, n);

        // If d == n, we've found a cycle but no factor
        // Try again with different parameters
        if &d == n {
            return None;
        }
    }

    Some(d)
}

/// Pollard's rho with Brent's improvement.
///
/// Brent's algorithm uses a different cycle detection that is about 24% faster.
/// It also batches GCD computations for additional speedup.
pub fn pollard_rho_brent(n: &Integer) -> Option<Integer> {
    if n.is_one() {
        return None;
    }

    // Check if n is even
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return Some(Integer::from(2i64));
    }

    // Check if n is a perfect square
    let sqrt_n = integer_sqrt(n);
    if &(sqrt_n.clone() * sqrt_n.clone()) == n {
        return Some(sqrt_n);
    }

    let mut rng = ChaCha8Rng::seed_from_u64(0xDEADBEEF);

    // Try with different parameters
    for _ in 0..20 {
        let c = Integer::from(rng.gen_range(1i64..1000));
        let y0 = Integer::from(rng.gen_range(2i64..1000));

        if let Some(factor) = pollard_rho_brent_with_params(n, y0, c) {
            if !factor.is_one() && &factor != n {
                return Some(factor);
            }
        }
    }

    None
}

/// Brent's rho algorithm with specific parameters.
fn pollard_rho_brent_with_params(n: &Integer, y0: Integer, c: Integer) -> Option<Integer> {
    let mut y = y0;
    let mut r = Integer::one();
    let mut q = Integer::one();
    let mut g = Integer::one();
    let mut ys = Integer::zero();
    let mut x = Integer::zero();

    let m = Integer::from(128i64); // Batch size for GCD

    let f = |x: Integer| -> Integer { ((x.clone() * x.clone()) + c.clone()) % n.clone() };

    while g.is_one() {
        x = y.clone();

        // Move y forward by r steps
        let r_val = r.to_i64().unwrap_or(1000) as usize;
        for _ in 0..r_val {
            y = f(y);
        }

        let mut k = Integer::zero();
        while &k < &r && g.is_one() {
            ys = y.clone();

            // Batch GCD computation
            let steps = std::cmp::min(
                m.to_i64().unwrap_or(128) as usize,
                (r.clone() - k.clone()).to_i64().unwrap_or(128) as usize,
            );

            for _ in 0..steps {
                y = f(y);
                let diff = abs_diff(&x, &y);
                q = (q.clone() * diff) % n.clone();
            }

            g = gcd(&q, n);
            k = k + Integer::from(steps as i64);
        }

        r = r.clone() * Integer::from(2i64);

        // Safety check: don't run forever
        if r.bit_len() > 64 {
            return None;
        }
    }

    // If g == n, we need to backtrack
    if &g == n {
        loop {
            ys = f(ys);
            g = gcd(&abs_diff(&x, &ys), n);
            if !g.is_one() {
                break;
            }
        }
    }

    if &g != n && !g.is_one() {
        Some(g)
    } else {
        None
    }
}

/// Compute |a - b|.
fn abs_diff(a: &Integer, b: &Integer) -> Integer {
    if a >= b {
        a.clone() - b.clone()
    } else {
        b.clone() - a.clone()
    }
}

/// Compute gcd(a, b) using Euclidean algorithm.
fn gcd(a: &Integer, b: &Integer) -> Integer {
    a.gcd(b)
}

/// Integer square root.
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
    fn test_pollard_rho_small() {
        // 15 = 3 × 5
        let factor = pollard_rho(&Integer::from(15i64));
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!(f == Integer::from(3i64) || f == Integer::from(5i64));
    }

    #[test]
    fn test_pollard_rho_semiprime() {
        // 1000003 × 1000033 = 1000036000099
        let n = Integer::from(1000003i64) * Integer::from(1000033i64);
        let factor = pollard_rho_brent(&n);
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!(f == Integer::from(1000003i64) || f == Integer::from(1000033i64));
    }

    #[test]
    fn test_pollard_rho_even() {
        let factor = pollard_rho(&Integer::from(100i64));
        assert_eq!(factor, Some(Integer::from(2i64)));
    }

    #[test]
    fn test_pollard_rho_brent() {
        // 1000000007 × 1000000009
        let n = Integer::from(1000000007i64) * Integer::from(1000000009i64);
        let factor = pollard_rho_brent(&n);
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!(f == Integer::from(1000000007i64) || f == Integer::from(1000000009i64));
    }

    #[test]
    fn test_pollard_rho_perfect_square() {
        let n = Integer::from(49i64); // 7^2
        let factor = pollard_rho_brent(&n);
        assert_eq!(factor, Some(Integer::from(7i64)));
    }
}
