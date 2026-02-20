//! Arithmetic functions.
//!
//! This module provides classical number-theoretic functions:
//! - Euler's totient φ(n)
//! - Carmichael's function λ(n)
//! - Divisor functions σ_k(n), τ(n)
//! - Möbius function μ(n)
//! - Radical rad(n)

use num_traits::{One, Zero};
use tertius_integers::Integer;

use crate::factorization::{factor, Factorization};

/// Euler's totient function φ(n).
///
/// Returns the count of integers in [1, n] that are coprime to n.
///
/// φ(n) = n × ∏(1 - 1/p) for all prime factors p of n
///
/// # Examples
///
/// ```
/// use tertius_numtheory::euler_phi;
/// use tertius_integers::Integer;
///
/// assert_eq!(euler_phi(&Integer::from(12i64)), Integer::from(4i64)); // {1, 5, 7, 11}
/// assert_eq!(euler_phi(&Integer::from(7i64)), Integer::from(6i64));  // 7 is prime
/// ```
pub fn euler_phi(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }
    if n.is_one() {
        return Integer::one();
    }

    let factorization = factor(&n.abs());
    euler_phi_from_factorization(&factorization)
}

/// Compute Euler's totient from a factorization.
pub fn euler_phi_from_factorization(f: &Factorization) -> Integer {
    let mut result = Integer::one();

    for (p, &e) in f.iter() {
        // φ(p^e) = p^(e-1) × (p - 1)
        let p_power = p.pow(e - 1);
        let p_minus_1 = p.clone() - Integer::one();
        result = result * p_power * p_minus_1;
    }

    result
}

/// Carmichael's function λ(n).
///
/// Returns the smallest positive integer m such that a^m ≡ 1 (mod n)
/// for all a coprime to n.
///
/// λ(n) = lcm(λ(p₁^e₁), λ(p₂^e₂), ...)
///
/// # Examples
///
/// ```
/// use tertius_numtheory::carmichael_lambda;
/// use tertius_integers::Integer;
///
/// assert_eq!(carmichael_lambda(&Integer::from(8i64)), Integer::from(2i64));
/// assert_eq!(carmichael_lambda(&Integer::from(15i64)), Integer::from(4i64));
/// ```
pub fn carmichael_lambda(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }
    if n.is_one() {
        return Integer::one();
    }

    let factorization = factor(&n.abs());
    carmichael_lambda_from_factorization(&factorization)
}

/// Compute Carmichael's function from a factorization.
pub fn carmichael_lambda_from_factorization(f: &Factorization) -> Integer {
    let mut result = Integer::one();

    for (p, &e) in f.iter() {
        let lambda_pe = if p == &Integer::from(2i64) && e >= 3 {
            // λ(2^e) = 2^(e-2) for e ≥ 3
            Integer::from(2i64).pow(e - 2)
        } else {
            // λ(p^e) = φ(p^e) = p^(e-1)(p-1) for odd p or 2^1, 2^2
            let p_power = p.pow(e - 1);
            let p_minus_1 = p.clone() - Integer::one();
            p_power * p_minus_1
        };

        result = lcm(&result, &lambda_pe);
    }

    result
}

/// List all divisors of n.
///
/// Returns divisors in sorted order.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::divisors;
/// use tertius_integers::Integer;
///
/// let divs = divisors(&Integer::from(12i64));
/// assert_eq!(divs, vec![1, 2, 3, 4, 6, 12].into_iter().map(Integer::from).collect::<Vec<_>>());
/// ```
pub fn divisors(n: &Integer) -> Vec<Integer> {
    if n.is_zero() {
        return vec![];
    }
    if n.is_one() {
        return vec![Integer::one()];
    }

    let factorization = factor(&n.abs());
    divisors_from_factorization(&factorization)
}

/// Compute divisors from a factorization.
pub fn divisors_from_factorization(f: &Factorization) -> Vec<Integer> {
    let mut result = vec![Integer::one()];

    for (p, &e) in f.iter() {
        let mut new_divisors = Vec::new();
        let mut power = Integer::one();

        for _ in 0..=e {
            for d in &result {
                new_divisors.push(d.clone() * power.clone());
            }
            power = power * p.clone();
        }

        result = new_divisors;
    }

    result.sort();
    result
}

/// Count of divisors τ(n) = σ₀(n).
///
/// # Examples
///
/// ```
/// use tertius_numtheory::divisor_count;
/// use tertius_integers::Integer;
///
/// assert_eq!(divisor_count(&Integer::from(12i64)), 6); // 1, 2, 3, 4, 6, 12
/// ```
pub fn divisor_count(n: &Integer) -> u64 {
    if n.is_zero() {
        return 0;
    }
    if n.is_one() {
        return 1;
    }

    let factorization = factor(&n.abs());
    divisor_count_from_factorization(&factorization)
}

/// Compute divisor count from a factorization.
pub fn divisor_count_from_factorization(f: &Factorization) -> u64 {
    f.iter().map(|(_, &e)| (e + 1) as u64).product()
}

/// Sum of divisors σ(n) = σ₁(n).
///
/// # Examples
///
/// ```
/// use tertius_numtheory::divisor_sum;
/// use tertius_integers::Integer;
///
/// assert_eq!(divisor_sum(&Integer::from(12i64)), Integer::from(28i64)); // 1+2+3+4+6+12
/// ```
pub fn divisor_sum(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }
    if n.is_one() {
        return Integer::one();
    }

    let factorization = factor(&n.abs());
    divisor_sum_from_factorization(&factorization)
}

/// Compute divisor sum from a factorization.
pub fn divisor_sum_from_factorization(f: &Factorization) -> Integer {
    let mut result = Integer::one();

    for (p, &e) in f.iter() {
        // σ(p^e) = (p^(e+1) - 1) / (p - 1)
        let p_power = p.pow(e + 1);
        let numerator = p_power - Integer::one();
        let denominator = p.clone() - Integer::one();
        result = result * (numerator / denominator);
    }

    result
}

/// Möbius function μ(n).
///
/// Returns:
/// - 1 if n is squarefree with an even number of prime factors
/// - -1 if n is squarefree with an odd number of prime factors
/// - 0 if n has a squared prime factor
///
/// # Examples
///
/// ```
/// use tertius_numtheory::mobius;
/// use tertius_integers::Integer;
///
/// assert_eq!(mobius(&Integer::from(1i64)), 1);
/// assert_eq!(mobius(&Integer::from(6i64)), 1);  // 2×3, even number of factors
/// assert_eq!(mobius(&Integer::from(30i64)), -1); // 2×3×5, odd number of factors
/// assert_eq!(mobius(&Integer::from(4i64)), 0);   // 2² has a square
/// ```
pub fn mobius(n: &Integer) -> i8 {
    if n.is_zero() {
        return 0;
    }
    if n.is_one() {
        return 1;
    }

    let factorization = factor(&n.abs());
    mobius_from_factorization(&factorization)
}

/// Compute Möbius function from a factorization.
pub fn mobius_from_factorization(f: &Factorization) -> i8 {
    // Check if squarefree
    for (_, &e) in f.iter() {
        if e > 1 {
            return 0;
        }
    }

    // (-1)^(number of prime factors)
    if f.num_distinct_factors() % 2 == 0 {
        1
    } else {
        -1
    }
}

/// Radical of n: rad(n) = ∏_{p|n} p.
///
/// The radical is the product of distinct prime factors.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::radical;
/// use tertius_integers::Integer;
///
/// assert_eq!(radical(&Integer::from(12i64)), Integer::from(6i64)); // 2×3
/// assert_eq!(radical(&Integer::from(72i64)), Integer::from(6i64)); // 72 = 2³×3², rad = 2×3
/// ```
pub fn radical(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }
    if n.is_one() {
        return Integer::one();
    }

    let factorization = factor(&n.abs());
    radical_from_factorization(&factorization)
}

/// Compute radical from a factorization.
pub fn radical_from_factorization(f: &Factorization) -> Integer {
    let mut result = Integer::one();
    for (p, _) in f.iter() {
        result = result * p.clone();
    }
    result
}

/// Compute lcm(a, b).
fn lcm(a: &Integer, b: &Integer) -> Integer {
    if a.is_zero() || b.is_zero() {
        return Integer::zero();
    }
    let g = a.gcd(b);
    (a.clone() / g) * b.clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_phi() {
        assert_eq!(euler_phi(&Integer::from(1i64)), Integer::from(1i64));
        assert_eq!(euler_phi(&Integer::from(2i64)), Integer::from(1i64));
        assert_eq!(euler_phi(&Integer::from(6i64)), Integer::from(2i64));
        assert_eq!(euler_phi(&Integer::from(12i64)), Integer::from(4i64));
        assert_eq!(euler_phi(&Integer::from(100i64)), Integer::from(40i64));
    }

    #[test]
    fn test_carmichael() {
        assert_eq!(carmichael_lambda(&Integer::from(1i64)), Integer::from(1i64));
        assert_eq!(carmichael_lambda(&Integer::from(8i64)), Integer::from(2i64));
        assert_eq!(
            carmichael_lambda(&Integer::from(15i64)),
            Integer::from(4i64)
        );
    }

    #[test]
    fn test_divisor_count() {
        assert_eq!(divisor_count(&Integer::from(1i64)), 1);
        assert_eq!(divisor_count(&Integer::from(12i64)), 6);
        assert_eq!(divisor_count(&Integer::from(100i64)), 9);
    }

    #[test]
    fn test_divisor_sum() {
        assert_eq!(divisor_sum(&Integer::from(1i64)), Integer::from(1i64));
        assert_eq!(divisor_sum(&Integer::from(12i64)), Integer::from(28i64));
    }

    #[test]
    fn test_mobius() {
        assert_eq!(mobius(&Integer::from(1i64)), 1);
        assert_eq!(mobius(&Integer::from(2i64)), -1);
        assert_eq!(mobius(&Integer::from(6i64)), 1);
        assert_eq!(mobius(&Integer::from(4i64)), 0);
        assert_eq!(mobius(&Integer::from(30i64)), -1);
    }

    #[test]
    fn test_radical() {
        assert_eq!(radical(&Integer::from(1i64)), Integer::from(1i64));
        assert_eq!(radical(&Integer::from(12i64)), Integer::from(6i64));
        assert_eq!(radical(&Integer::from(72i64)), Integer::from(6i64));
    }
}
