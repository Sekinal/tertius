//! Modular arithmetic and congruences.
//!
//! This module provides:
//! - Chinese Remainder Theorem
//! - Legendre and Jacobi symbols
//! - Tonelli-Shanks algorithm for modular square roots

use num_traits::{One, Zero};
use tertius_integers::Integer;

/// Chinese Remainder Theorem.
///
/// Given systems of congruences x ≡ aᵢ (mod nᵢ) where the nᵢ are pairwise coprime,
/// finds the unique solution x modulo N = ∏nᵢ.
///
/// # Arguments
///
/// * `residues` - The residues aᵢ
/// * `moduli` - The moduli nᵢ (must be pairwise coprime)
///
/// # Returns
///
/// Returns `Some((x, N))` where x is the solution and N is the product of moduli,
/// or `None` if the moduli are not pairwise coprime.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::crt;
/// use tertius_integers::Integer;
///
/// // x ≡ 2 (mod 3), x ≡ 3 (mod 5) => x ≡ 8 (mod 15)
/// let result = crt(
///     &[Integer::from(2i64), Integer::from(3i64)],
///     &[Integer::from(3i64), Integer::from(5i64)]
/// );
/// assert_eq!(result, Some((Integer::from(8i64), Integer::from(15i64))));
/// ```
pub fn crt(residues: &[Integer], moduli: &[Integer]) -> Option<(Integer, Integer)> {
    if residues.len() != moduli.len() {
        return None;
    }
    if residues.is_empty() {
        return Some((Integer::zero(), Integer::one()));
    }

    let mut result = residues[0].clone();
    let mut mod_product = moduli[0].clone();

    for i in 1..residues.len() {
        let (g, u, _) = extended_gcd(&mod_product, &moduli[i]);

        if !g.is_one() {
            // Check if the system is solvable
            let diff = residues[i].clone() - result.clone();
            if !(diff.clone() % g.clone()).is_zero() {
                return None;
            }
        }

        // result += u * (residues[i] - result) * (mod_product / g)
        let diff = residues[i].clone() - result.clone();
        let adjustment = (u * diff) % moduli[i].clone();
        result = result + adjustment * mod_product.clone();

        mod_product = mod_product * moduli[i].clone() / g;
        result = ((result % mod_product.clone()) + mod_product.clone()) % mod_product.clone();
    }

    Some((result, mod_product))
}

/// Legendre symbol (a/p) for prime p.
///
/// Returns:
/// - 1 if a is a quadratic residue mod p
/// - -1 if a is a non-residue mod p
/// - 0 if p divides a
///
/// # Examples
///
/// ```
/// use tertius_numtheory::legendre_symbol;
/// use tertius_integers::Integer;
///
/// assert_eq!(legendre_symbol(&Integer::from(2i64), &Integer::from(7i64)), 1);
/// assert_eq!(legendre_symbol(&Integer::from(3i64), &Integer::from(7i64)), -1);
/// ```
pub fn legendre_symbol(a: &Integer, p: &Integer) -> i32 {
    // Use Euler's criterion: (a/p) ≡ a^((p-1)/2) (mod p)
    let a_mod_p = ((a.clone() % p.clone()) + p.clone()) % p.clone();

    if a_mod_p.is_zero() {
        return 0;
    }

    let exp = (p.clone() - Integer::one()) / Integer::from(2i64);
    let result = mod_pow(&a_mod_p, &exp, p);

    if result.is_one() {
        1
    } else if result == p.clone() - Integer::one() {
        -1
    } else {
        0
    }
}

/// Jacobi symbol (a/n) for odd n.
///
/// Generalization of Legendre symbol to composite odd n.
/// Computed using quadratic reciprocity.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::jacobi_symbol;
/// use tertius_integers::Integer;
///
/// assert_eq!(jacobi_symbol(&Integer::from(5i64), &Integer::from(21i64)), 1);
/// ```
pub fn jacobi_symbol(a: &Integer, n: &Integer) -> i32 {
    if n.is_one() {
        return 1;
    }

    let mut a = ((a.clone() % n.clone()) + n.clone()) % n.clone();
    let mut n = n.clone();
    let mut result = 1i32;

    while !a.is_zero() {
        // Factor out powers of 2
        let mut twos = 0u32;
        while (a.clone() % Integer::from(2i64)).is_zero() {
            a = a / Integer::from(2i64);
            twos += 1;
        }

        // Apply quadratic reciprocity for 2
        if twos % 2 == 1 {
            let n_mod_8 = (n.clone() % Integer::from(8i64)).to_i64().unwrap_or(0);
            if n_mod_8 == 3 || n_mod_8 == 5 {
                result = -result;
            }
        }

        // Swap a and n
        std::mem::swap(&mut a, &mut n);

        // Apply quadratic reciprocity
        let a_mod_4 = (a.clone() % Integer::from(4i64)).to_i64().unwrap_or(0);
        let n_mod_4 = (n.clone() % Integer::from(4i64)).to_i64().unwrap_or(0);
        if a_mod_4 == 3 && n_mod_4 == 3 {
            result = -result;
        }

        a = a % n.clone();
    }

    if n.is_one() {
        result
    } else {
        0
    }
}

/// Tonelli-Shanks algorithm for modular square roots.
///
/// Finds x such that x² ≡ a (mod p) where p is an odd prime.
///
/// # Returns
///
/// Returns `Some(x)` if a is a quadratic residue mod p, `None` otherwise.
/// Note: if x is a solution, so is p - x.
///
/// # Examples
///
/// ```
/// use tertius_numtheory::tonelli_shanks;
/// use tertius_integers::Integer;
///
/// let x = tonelli_shanks(&Integer::from(2i64), &Integer::from(7i64)).unwrap();
/// assert!((x.clone() * x.clone()) % Integer::from(7i64) == Integer::from(2i64));
/// ```
pub fn tonelli_shanks(a: &Integer, p: &Integer) -> Option<Integer> {
    let a = ((a.clone() % p.clone()) + p.clone()) % p.clone();

    if a.is_zero() {
        return Some(Integer::zero());
    }

    // Check if a is a quadratic residue
    if legendre_symbol(&a, p) != 1 {
        return None;
    }

    // Special case: p ≡ 3 (mod 4)
    let p_mod_4 = (p.clone() % Integer::from(4i64)).to_i64().unwrap_or(0);
    if p_mod_4 == 3 {
        let exp = (p.clone() + Integer::one()) / Integer::from(4i64);
        return Some(mod_pow(&a, &exp, p));
    }

    // Factor p - 1 = Q * 2^S where Q is odd
    let mut q = p.clone() - Integer::one();
    let mut s = 0u32;
    while (q.clone() % Integer::from(2i64)).is_zero() {
        q = q / Integer::from(2i64);
        s += 1;
    }

    // Find a quadratic non-residue z
    let mut z = Integer::from(2i64);
    while legendre_symbol(&z, p) != -1 {
        z = z + Integer::one();
    }

    let mut m = s;
    let mut c = mod_pow(&z, &q, p);
    let mut t = mod_pow(&a, &q, p);
    let mut r = mod_pow(&a, &((q.clone() + Integer::one()) / Integer::from(2i64)), p);

    loop {
        if t.is_one() {
            return Some(r);
        }

        // Find least i such that t^(2^i) = 1
        let mut i = 1u32;
        let mut temp = (t.clone() * t.clone()) % p.clone();
        while !temp.is_one() {
            temp = (temp.clone() * temp.clone()) % p.clone();
            i += 1;
            if i >= m {
                return None; // Should not happen if a is a QR
            }
        }

        // b = c^(2^(m-i-1))
        let exp = Integer::from(2i64).pow(m - i - 1);
        let b = mod_pow(&c, &exp, p);

        m = i;
        c = (b.clone() * b.clone()) % p.clone();
        t = (t * c.clone()) % p.clone();
        r = (r * b) % p.clone();
    }
}

/// Extended Euclidean algorithm.
///
/// Returns (gcd, x, y) such that gcd = a*x + b*y.
pub fn extended_gcd(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    if b.is_zero() {
        return (a.clone(), Integer::one(), Integer::zero());
    }

    let (g, x1, y1) = extended_gcd(b, &(a.clone() % b.clone()));
    let x = y1.clone();
    let y = x1 - (a.clone() / b.clone()) * y1;
    (g, x, y)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crt_simple() {
        // x ≡ 2 (mod 3), x ≡ 3 (mod 5) => x ≡ 8 (mod 15)
        let result = crt(
            &[Integer::from(2i64), Integer::from(3i64)],
            &[Integer::from(3i64), Integer::from(5i64)],
        );
        let (x, m) = result.unwrap();
        assert_eq!(m, Integer::from(15i64));
        assert_eq!(x.clone() % Integer::from(3i64), Integer::from(2i64));
        assert_eq!(x % Integer::from(5i64), Integer::from(3i64));
    }

    #[test]
    fn test_legendre() {
        // 2 is a QR mod 7 (3² = 9 ≡ 2)
        assert_eq!(
            legendre_symbol(&Integer::from(2i64), &Integer::from(7i64)),
            1
        );
        // 3 is not a QR mod 7
        assert_eq!(
            legendre_symbol(&Integer::from(3i64), &Integer::from(7i64)),
            -1
        );
    }

    #[test]
    fn test_jacobi() {
        // (5/21) = (5/3)(5/7) = (-1)(-1) = 1
        assert_eq!(
            jacobi_symbol(&Integer::from(5i64), &Integer::from(21i64)),
            1
        );
    }

    #[test]
    fn test_tonelli_shanks() {
        // √2 mod 7 = 3 or 4 (since 3² = 9 ≡ 2, 4² = 16 ≡ 2)
        let x = tonelli_shanks(&Integer::from(2i64), &Integer::from(7i64)).unwrap();
        assert!((x.clone() * x) % Integer::from(7i64) == Integer::from(2i64));
    }

    #[test]
    fn test_tonelli_shanks_no_solution() {
        // 3 is not a QR mod 7
        assert!(tonelli_shanks(&Integer::from(3i64), &Integer::from(7i64)).is_none());
    }
}
