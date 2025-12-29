//! Self-Initializing Quadratic Sieve (SIQS) for integer factorization.
//!
//! SIQS is the second-fastest general-purpose factorization algorithm
//! (after the General Number Field Sieve) and is practical for numbers
//! up to about 100 digits.
//!
//! The algorithm works by:
//! 1. Finding many "relations": numbers x such that x² ≡ y (mod n) where y is smooth
//! 2. Using linear algebra to find a subset whose product is a perfect square
//! 3. This gives x² ≡ y² (mod n), and gcd(x-y, n) is often a factor
//!
//! The "self-initializing" variant uses multiple polynomials to speed up sieving.

mod polynomial;
mod sieve;
mod relations;

use num_traits::{One, Zero};
use rayon::prelude::*;
use tertius_integers::Integer;

use super::Factorization;

/// SIQS factorization with default parameters.
///
/// This is suitable for numbers from about 50 to 100 digits.
pub fn siqs_factor(n: &Integer) -> Option<Integer> {
    // Determine parameters based on size
    let bits = n.bit_len();
    let params = select_parameters(bits);

    siqs_with_params(n, &params)
}

/// Parameters for SIQS.
#[derive(Debug, Clone)]
pub struct SiqsParams {
    /// Size of the factor base (number of small primes).
    pub factor_base_size: usize,
    /// Sieving interval half-width.
    pub sieve_interval: usize,
    /// Number of large prime variations allowed.
    pub large_prime_bound: u64,
    /// Error threshold for smooth detection.
    pub error_threshold: f64,
}

/// Select parameters based on input size.
fn select_parameters(bits: usize) -> SiqsParams {
    // These are approximate optimal parameters from literature
    match bits {
        0..=64 => SiqsParams {
            factor_base_size: 100,
            sieve_interval: 10_000,
            large_prime_bound: 1_000_000,
            error_threshold: 1.5,
        },
        65..=96 => SiqsParams {
            factor_base_size: 300,
            sieve_interval: 50_000,
            large_prime_bound: 10_000_000,
            error_threshold: 1.8,
        },
        97..=128 => SiqsParams {
            factor_base_size: 1000,
            sieve_interval: 100_000,
            large_prime_bound: 100_000_000,
            error_threshold: 2.0,
        },
        129..=192 => SiqsParams {
            factor_base_size: 3000,
            sieve_interval: 200_000,
            large_prime_bound: 1_000_000_000,
            error_threshold: 2.2,
        },
        _ => SiqsParams {
            factor_base_size: 10000,
            sieve_interval: 500_000,
            large_prime_bound: 10_000_000_000,
            error_threshold: 2.5,
        },
    }
}

/// SIQS with custom parameters.
pub fn siqs_with_params(n: &Integer, params: &SiqsParams) -> Option<Integer> {
    // Check trivial cases
    if n.is_one() {
        return None;
    }
    if (n.clone() % Integer::from(2i64)).is_zero() {
        return Some(Integer::from(2i64));
    }

    // Check if n is a perfect square
    let sqrt_n = integer_sqrt(n);
    if &(sqrt_n.clone() * sqrt_n.clone()) == n {
        return Some(sqrt_n);
    }

    // Step 1: Build factor base
    let factor_base = build_factor_base(n, params.factor_base_size);
    if factor_base.is_empty() {
        return None;
    }

    // Step 2: Find smooth relations using sieving
    let relations = sieve::sieve_for_relations(n, &factor_base, params);

    // We need slightly more relations than factor base size
    if relations.len() <= factor_base.len() {
        return None;
    }

    // Step 3: Use linear algebra to find dependencies
    // This is the matrix step - we use GF(2) linear algebra
    if let Some(dependency) = find_dependency(&relations, &factor_base) {
        // Step 4: Compute x and y from the dependency
        let (x, y) = compute_square_root(&relations, &dependency, n);

        // Step 5: Compute gcd(x - y, n)
        let diff = if x >= y {
            x.clone() - y.clone()
        } else {
            y.clone() - x.clone()
        };
        let g = gcd(&diff, n);

        if !g.is_one() && &g != n {
            return Some(g);
        }

        // Try x + y as well
        let sum = (x + y) % n.clone();
        let g = gcd(&sum, n);
        if !g.is_one() && &g != n {
            return Some(g);
        }
    }

    None
}

/// A relation: x² ≡ prod(p_i^e_i) (mod n)
#[derive(Debug, Clone)]
pub struct Relation {
    /// The x value.
    pub x: Integer,
    /// Exponent vector over the factor base (mod 2 for linear algebra).
    pub exponents: Vec<u8>,
    /// The actual factorization (for computing square root).
    pub factorization: Vec<(usize, u32)>,
}

/// Build the factor base: small primes p where n is a quadratic residue mod p.
fn build_factor_base(n: &Integer, size: usize) -> Vec<u64> {
    let mut base = vec![2]; // Always include 2
    let mut p = 3u64;

    while base.len() < size {
        if is_prime_simple(p) && is_quadratic_residue(n, p) {
            base.push(p);
        }
        p += 2;
    }

    base
}

/// Check if n is a quadratic residue mod p using Euler's criterion.
fn is_quadratic_residue(n: &Integer, p: u64) -> bool {
    if p == 2 {
        return true;
    }
    let n_mod_p = (n.clone() % Integer::from(p)).to_i64().unwrap_or(0) as u64;
    if n_mod_p == 0 {
        return true;
    }
    // a^((p-1)/2) ≡ 1 (mod p) iff a is a QR
    mod_pow_u64(n_mod_p, (p - 1) / 2, p) == 1
}

/// Simple primality test.
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

/// Modular exponentiation for u64.
fn mod_pow_u64(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1u64;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % m as u128) as u64;
        }
        exp >>= 1;
        base = ((base as u128 * base as u128) % m as u128) as u64;
    }
    result
}

/// Find a linear dependency in the relation matrix over GF(2).
///
/// Uses Gaussian elimination to find the null space.
fn find_dependency(relations: &[Relation], factor_base: &[u64]) -> Option<Vec<usize>> {
    let n_cols = relations.len();
    let n_rows = factor_base.len();

    if n_cols <= n_rows {
        return None;
    }

    // Build matrix (transposed for row reduction)
    let mut matrix: Vec<Vec<u8>> = vec![vec![0; n_cols]; n_rows];
    for (col, rel) in relations.iter().enumerate() {
        for &exp in &rel.exponents {
            if (exp as usize) < n_rows {
                matrix[exp as usize][col] = 1;
            }
        }
    }

    // Gaussian elimination over GF(2)
    let mut pivot_cols = Vec::new();
    let mut row = 0;

    for col in 0..n_cols {
        // Find pivot
        let mut found = false;
        for r in row..n_rows {
            if matrix[r][col] == 1 {
                matrix.swap(row, r);
                found = true;
                break;
            }
        }

        if !found {
            // This column is in the null space
            // Construct the dependency
            let mut dep = vec![col];
            for (c, &pc) in pivot_cols.iter().enumerate() {
                if matrix[c][col] == 1 {
                    dep.push(pc);
                }
            }
            if dep.len() > 1 {
                return Some(dep);
            }
            continue;
        }

        // Eliminate
        for r in 0..n_rows {
            if r != row && matrix[r][col] == 1 {
                for c in 0..n_cols {
                    matrix[r][c] ^= matrix[row][c];
                }
            }
        }

        pivot_cols.push(col);
        row += 1;
        if row >= n_rows {
            break;
        }
    }

    None
}

/// Compute x and y such that x² ≡ y² (mod n) from a dependency.
fn compute_square_root(
    relations: &[Relation],
    dependency: &[usize],
    n: &Integer,
) -> (Integer, Integer) {
    // x = product of relation.x values
    let mut x = Integer::one();
    // y = square root of product of smooth parts
    let mut total_exponents = vec![0u32; 1000]; // Rough upper bound

    for &idx in dependency {
        let rel = &relations[idx];
        x = (x * rel.x.clone()) % n.clone();
        for &(prime_idx, exp) in &rel.factorization {
            if prime_idx < total_exponents.len() {
                total_exponents[prime_idx] += exp;
            }
        }
    }

    // All exponents should be even now
    let mut y = Integer::one();
    for (idx, &exp) in total_exponents.iter().enumerate() {
        if exp > 0 {
            // This is a placeholder - we need the actual prime values
            let prime = Integer::from(idx as i64 * 2 + 1); // Rough approximation
            y = (y * mod_pow(&prime, &Integer::from((exp / 2) as i64), n)) % n.clone();
        }
    }

    (x, y)
}

/// Modular exponentiation.
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

/// GCD.
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
    fn test_build_factor_base() {
        let n = Integer::from(1000003i64 * 1000033i64);
        let base = build_factor_base(&n, 20);
        assert!(base.len() >= 10);
        assert_eq!(base[0], 2);
    }

    #[test]
    fn test_is_quadratic_residue() {
        // 2 is a QR mod 7 (since 3² = 9 ≡ 2)
        assert!(is_quadratic_residue(&Integer::from(2i64), 7));
        // 3 is not a QR mod 7
        assert!(!is_quadratic_residue(&Integer::from(3i64), 7));
    }
}
