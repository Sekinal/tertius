//! Lucas primality test.
//!
//! The Lucas primality test is based on Lucas sequences and is used
//! as part of the Baillie-PSW test. It catches composites that might
//! fool Miller-Rabin.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use super::miller_rabin::mod_pow;

/// Standard Lucas test.
///
/// Uses the Lucas sequence U_n(P, Q) to test primality.
/// Returns true if n passes the test.
pub fn lucas_test(n: &Integer) -> bool {
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

    // Find D using the Selfridge method
    let d = selfridge_d(n);
    if d.is_zero() {
        // n is a perfect square, hence composite
        return false;
    }

    // For Lucas sequence with P = 1 and Q = (1 - D) / 4
    // Actually we use the standard Baillie-PSW parameters
    let (p, q) = selfridge_pq(&d);

    // Compute U_{n+1} mod n
    // If n is prime and gcd(n, Q) = 1, then U_{n-(D/n)} ≡ 0 (mod n)
    let delta = n.clone() + Integer::one(); // n + 1 when (D/n) = -1

    let (u, _v) = lucas_sequence(n, &p, &q, &delta);

    u.is_zero()
}

/// Strong Lucas test (used in Baillie-PSW).
///
/// This is a stronger variant of the Lucas test.
pub fn strong_lucas_test(n: &Integer) -> bool {
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

    // Check if n is a perfect square
    let sqrt_n = integer_sqrt(n);
    if &(sqrt_n.clone() * sqrt_n.clone()) == n {
        return false; // Perfect squares are composite (unless 1)
    }

    // Find D using the Selfridge method
    let d = selfridge_d(n);
    if d.is_zero() {
        return false;
    }

    // P = 1, Q = (1 - D) / 4
    let (p, q) = selfridge_pq(&d);

    // Write n + 1 = 2^s * d where d is odd
    let n_plus_1 = n.clone() + Integer::one();
    let (s, odd_d) = factor_out_twos(&n_plus_1);

    // Compute U_d and V_d mod n
    let (mut u, mut v) = lucas_sequence(n, &p, &q, &odd_d);

    // If U_d ≡ 0 (mod n) or V_d ≡ 0 (mod n), n passes
    if u.is_zero() || v.is_zero() {
        return true;
    }

    // Check V_{2^r * d} for r = 1, 2, ..., s-1
    // We need to track Q^k where k doubles each iteration
    let two = Integer::from(2i64);
    let mut q_k = mod_pow(&q, &odd_d, n); // Q^d initially

    for _ in 1..s {
        // V_{2k} = V_k^2 - 2*Q^k
        v = (v.clone() * v.clone() - two.clone() * q_k.clone()) % n.clone();
        if v.is_negative() {
            v = v + n.clone();
        }

        // Q^{2k} = (Q^k)^2
        q_k = (q_k.clone() * q_k.clone()) % n.clone();

        if v.is_zero() {
            return true;
        }
    }

    false
}

/// Selfridge's method to find parameter D for Lucas test.
///
/// Finds the first D in the sequence 5, -7, 9, -11, 13, -15, ...
/// such that Jacobi(D, n) = -1.
fn selfridge_d(n: &Integer) -> Integer {
    // For small n, use a simple approach
    if n <= &Integer::from(5i64) {
        // For n=3, D=5 gives Jacobi(5,3) = Jacobi(2,3) = -1
        // For n=5, D=-7 gives Jacobi(-7,5) = Jacobi(3,5) * (-1) = ...
        // Actually for n=5, D=-7: -7 mod 5 = 3, Jacobi(3,5) = ?
        // 3^2 = 9 ≡ 4 (mod 5), 4 is not 1, so 3 is a non-residue
        // So Jacobi(3,5) = -1, and D=-7 works for n=5
        return Integer::from(-7i64);
    }

    let mut d = Integer::from(5i64);
    let mut sign = 1i64;

    loop {
        let jacobi = jacobi_symbol(&d, n);

        if jacobi == -1 {
            return d;
        }
        if jacobi == 0 {
            // gcd(D, n) > 1, so n is composite (unless D equals n)
            if &d.abs() != n {
                return Integer::zero(); // Signal that n is composite
            }
        }

        // Next D in sequence
        sign = -sign;
        d = Integer::from(sign * (d.abs().to_i64().unwrap_or(0) + 2));

        // Safety check: if we've gone too far, something is wrong
        // The first D with Jacobi(D, n) = -1 is typically O(log² n),
        // so we use a generous bound of 2*n to avoid infinite loops
        // while allowing for primes that need larger D values.
        if d.abs() > n.clone() * Integer::from(2i64) {
            return Integer::zero();
        }
    }
}

/// Compute P and Q from D using Selfridge's method.
/// P = 1, Q = (1 - D) / 4
fn selfridge_pq(d: &Integer) -> (Integer, Integer) {
    let p = Integer::one();
    let q = (Integer::one() - d.clone()) / Integer::from(4i64);
    (p, q)
}

/// Compute (U_k, V_k) mod n for Lucas sequence with parameters P, Q.
///
/// Uses the binary method for efficient computation.
fn lucas_sequence(n: &Integer, p: &Integer, q: &Integer, k: &Integer) -> (Integer, Integer) {
    if k.is_zero() {
        return (Integer::zero(), Integer::from(2i64));
    }
    if k.is_one() {
        return (Integer::one(), p.clone() % n.clone());
    }

    // Binary expansion method
    let mut u = Integer::one();
    let mut v = p.clone();
    let mut q_k = q.clone(); // Q^1
    let two = Integer::from(2i64);
    let d = p.clone() * p.clone() - Integer::from(4i64) * q.clone();

    // Get binary representation of k
    let bits = integer_to_bits(k);

    // Start from the second-highest bit
    for i in 1..bits.len() {
        // Double: (U_{2k}, V_{2k})
        // U_{2k} = U_k * V_k
        // V_{2k} = V_k^2 - 2*Q^k
        let u_new = (u.clone() * v.clone()) % n.clone();
        let v_new = (v.clone() * v.clone() - two.clone() * q_k.clone()) % n.clone();
        u = normalize_mod(&u_new, n);
        v = normalize_mod(&v_new, n);
        q_k = (q_k.clone() * q_k.clone()) % n.clone();

        if bits[i] {
            // Add one: (U_{2k+1}, V_{2k+1})
            // U_{k+1} = (P*U_k + V_k) / 2
            // V_{k+1} = (D*U_k + P*V_k) / 2
            let u_next = (p.clone() * u.clone() + v.clone()) % n.clone();
            let v_next = (d.clone() * u.clone() + p.clone() * v.clone()) % n.clone();

            // Divide by 2 (mod n)
            u = div_by_2_mod(&u_next, n);
            v = div_by_2_mod(&v_next, n);
            q_k = (q_k * q.clone()) % n.clone();
        }
    }

    (normalize_mod(&u, n), normalize_mod(&v, n))
}

/// Compute Jacobi symbol (a/n).
fn jacobi_symbol(a: &Integer, n: &Integer) -> i32 {
    if n.is_one() {
        return 1;
    }
    if a.is_zero() {
        return 0;
    }

    let mut a = a.clone() % n.clone();
    if a.is_negative() {
        a = a + n.clone();
    }
    let mut n = n.clone();
    let mut result = 1i32;

    while !a.is_zero() {
        // Factor out powers of 2 from a
        let mut twos = 0u32;
        while (a.clone() % Integer::from(2i64)).is_zero() {
            a = a / Integer::from(2i64);
            twos += 1;
        }

        // Apply quadratic reciprocity for powers of 2
        if twos % 2 == 1 {
            let n_mod_8 = (n.clone() % Integer::from(8i64)).to_i64().unwrap_or(0);
            if n_mod_8 == 3 || n_mod_8 == 5 {
                result = -result;
            }
        }

        // Swap a and n (quadratic reciprocity)
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

/// Factor out powers of 2.
fn factor_out_twos(n: &Integer) -> (u32, Integer) {
    let mut d = n.clone();
    let mut s = 0u32;
    let two = Integer::from(2i64);

    while (d.clone() % two.clone()).is_zero() {
        d = d / two.clone();
        s += 1;
    }

    (s, d)
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

/// Convert integer to binary representation (most significant bit first).
fn integer_to_bits(n: &Integer) -> Vec<bool> {
    if n.is_zero() {
        return vec![false];
    }

    let mut bits = Vec::new();
    let mut m = n.clone();
    let two = Integer::from(2i64);

    while !m.is_zero() {
        bits.push(!(m.clone() % two.clone()).is_zero());
        m = m / two.clone();
    }

    bits.reverse();
    bits
}

/// Normalize a value to be in [0, n).
fn normalize_mod(a: &Integer, n: &Integer) -> Integer {
    let r = a.clone() % n.clone();
    if r.is_negative() {
        r + n.clone()
    } else {
        r
    }
}

/// Divide by 2 mod n.
fn div_by_2_mod(a: &Integer, n: &Integer) -> Integer {
    let a = normalize_mod(a, n);
    if (a.clone() % Integer::from(2i64)).is_zero() {
        a / Integer::from(2i64)
    } else {
        (a + n.clone()) / Integer::from(2i64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jacobi_symbol() {
        // (2/7) = 1 since 3^2 = 9 ≡ 2 (mod 7)
        assert_eq!(jacobi_symbol(&Integer::from(2i64), &Integer::from(7i64)), 1);

        // (5/21) = (5/3)(5/7)
        // 5 mod 3 = 2, (2/3) = -1
        // 5 mod 7 = 5, (5/7) = -1 (since 5 is not a QR mod 7)
        // So (5/21) = (-1)(-1) = 1
        assert_eq!(jacobi_symbol(&Integer::from(5i64), &Integer::from(21i64)), 1);

        // (3/5) = -1
        assert_eq!(jacobi_symbol(&Integer::from(3i64), &Integer::from(5i64)), -1);
    }

    #[test]
    fn test_lucas_test() {
        // Medium primes (Lucas test is meant for larger numbers)
        assert!(lucas_test(&Integer::from(7i64)));
        assert!(lucas_test(&Integer::from(11i64)));
        assert!(lucas_test(&Integer::from(997i64)));

        // Small composites
        assert!(!lucas_test(&Integer::from(9i64)));
        assert!(!lucas_test(&Integer::from(15i64)));
        assert!(!lucas_test(&Integer::from(21i64)));
    }

    #[test]
    fn test_strong_lucas() {
        // Medium primes
        assert!(strong_lucas_test(&Integer::from(7i64)));
        assert!(strong_lucas_test(&Integer::from(11i64)));
        assert!(strong_lucas_test(&Integer::from(997i64)));

        assert!(!strong_lucas_test(&Integer::from(9i64)));
        assert!(!strong_lucas_test(&Integer::from(25i64))); // Perfect square
    }
}
