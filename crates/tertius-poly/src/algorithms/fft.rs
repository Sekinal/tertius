//! FFT-based polynomial multiplication using CRT.
//!
//! For integer polynomials, we use the Chinese Remainder Theorem (CRT)
//! with multiple NTT primes to reconstruct the exact result.

use num_traits::Zero;
use tertius_integers::{Integer, ModInt};

use crate::algorithms::ntt::primes;

/// Type aliases for the three NTT fields.
type F1 = ModInt<{ primes::P1 }>;
type F2 = ModInt<{ primes::P2 }>;
type F3 = ModInt<{ primes::P3 }>;

/// Computes the NTT for a given prime and root.
fn ntt_with_prime<const P: u64, const G: u64>(a: &mut [ModInt<P>]) {
    let n = a.len();
    if n == 1 {
        return;
    }

    // Bit-reversal permutation
    let log_n = n.trailing_zeros();
    for i in 0..n {
        let j = (i as u32).reverse_bits() >> (32 - log_n);
        if i < j as usize {
            a.swap(i, j as usize);
        }
    }

    // Cooley-Tukey iterative NTT
    let mut len = 2;
    while len <= n {
        let w_len = ModInt::<P>::new(G).pow((P - 1) / len as u64);
        for i in (0..n).step_by(len) {
            let mut w = ModInt::<P>::new(1);
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w = w * w_len;
            }
        }
        len *= 2;
    }
}

/// Computes the inverse NTT for a given prime and root.
fn intt_with_prime<const P: u64, const G: u64>(a: &mut [ModInt<P>]) {
    let n = a.len();
    if n == 1 {
        return;
    }

    // Bit-reversal permutation
    let log_n = n.trailing_zeros();
    for i in 0..n {
        let j = (i as u32).reverse_bits() >> (32 - log_n);
        if i < j as usize {
            a.swap(i, j as usize);
        }
    }

    // Cooley-Tukey iterative INTT
    let mut len = 2;
    while len <= n {
        let w_len = ModInt::<P>::new(G).pow((P - 1) / len as u64).inv().unwrap();
        for i in (0..n).step_by(len) {
            let mut w = ModInt::<P>::new(1);
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w = w * w_len;
            }
        }
        len *= 2;
    }

    // Divide by n
    let n_inv = ModInt::<P>::new(n as u64).inv().unwrap();
    for x in a.iter_mut() {
        *x = *x * n_inv;
    }
}

/// Multiplies two integer polynomials using multi-modular NTT and CRT.
///
/// This method:
/// 1. Converts coefficients to three NTT-friendly prime fields
/// 2. Performs NTT multiplication in each field
/// 3. Uses CRT to reconstruct the integer result
///
/// The result is exact for coefficients up to about 2^90.
pub fn fft_multiply_integers(a: &[Integer], b: &[Integer]) -> Vec<Integer> {
    if a.is_empty() || b.is_empty() {
        return vec![Integer::new(0)];
    }

    let result_len = a.len() + b.len() - 1;
    let n = result_len.next_power_of_two();

    // Convert to three prime fields
    let mut a1: Vec<F1> = a.iter().map(|x| to_mod::<{ primes::P1 }>(x)).collect();
    let mut a2: Vec<F2> = a.iter().map(|x| to_mod::<{ primes::P2 }>(x)).collect();
    let mut a3: Vec<F3> = a.iter().map(|x| to_mod::<{ primes::P3 }>(x)).collect();

    let mut b1: Vec<F1> = b.iter().map(|x| to_mod::<{ primes::P1 }>(x)).collect();
    let mut b2: Vec<F2> = b.iter().map(|x| to_mod::<{ primes::P2 }>(x)).collect();
    let mut b3: Vec<F3> = b.iter().map(|x| to_mod::<{ primes::P3 }>(x)).collect();

    // Pad to power of 2
    a1.resize(n, F1::new(0));
    a2.resize(n, F2::new(0));
    a3.resize(n, F3::new(0));
    b1.resize(n, F1::new(0));
    b2.resize(n, F2::new(0));
    b3.resize(n, F3::new(0));

    // Forward NTT
    ntt_with_prime::<{ primes::P1 }, { primes::G1 }>(&mut a1);
    ntt_with_prime::<{ primes::P2 }, { primes::G2 }>(&mut a2);
    ntt_with_prime::<{ primes::P3 }, { primes::G3 }>(&mut a3);
    ntt_with_prime::<{ primes::P1 }, { primes::G1 }>(&mut b1);
    ntt_with_prime::<{ primes::P2 }, { primes::G2 }>(&mut b2);
    ntt_with_prime::<{ primes::P3 }, { primes::G3 }>(&mut b3);

    // Pointwise multiplication
    for i in 0..n {
        a1[i] = a1[i] * b1[i];
        a2[i] = a2[i] * b2[i];
        a3[i] = a3[i] * b3[i];
    }

    // Inverse NTT
    intt_with_prime::<{ primes::P1 }, { primes::G1 }>(&mut a1);
    intt_with_prime::<{ primes::P2 }, { primes::G2 }>(&mut a2);
    intt_with_prime::<{ primes::P3 }, { primes::G3 }>(&mut a3);

    // CRT reconstruction
    let mut result = Vec::with_capacity(result_len);
    for i in 0..result_len {
        let r1 = a1[i].value();
        let r2 = a2[i].value();
        let r3 = a3[i].value();
        result.push(crt3(r1, r2, r3));
    }

    result
}

/// Converts an Integer to a modular integer.
fn to_mod<const P: u64>(x: &Integer) -> ModInt<P> {
    if let Some(v) = x.to_i64() {
        ModInt::from_signed(v)
    } else {
        // For large integers, compute x mod P
        let p = Integer::from(P as i64);
        let r = x.clone() % p;
        ModInt::from_signed(r.to_i64().unwrap_or(0))
    }
}

/// Reconstructs an integer from residues modulo three primes using CRT.
fn crt3(r1: u64, r2: u64, r3: u64) -> Integer {
    // M = P1 * P2 * P3
    // We compute: x = r1*M1*y1 + r2*M2*y2 + r3*M3*y3 (mod M)
    // where Mi = M/Pi and yi = Mi^(-1) (mod Pi)

    let p1 = Integer::from(primes::P1 as i64);
    let p2 = Integer::from(primes::P2 as i64);
    let p3 = Integer::from(primes::P3 as i64);

    let m12 = p1.clone() * p2.clone();
    let m = m12.clone() * p3.clone();

    // M1 = P2 * P3, y1 = M1^(-1) mod P1
    let m1 = p2.clone() * p3.clone();
    let y1 = mod_inverse(m1.clone() % p1.clone(), p1.clone());

    // M2 = P1 * P3, y2 = M2^(-1) mod P2
    let m2 = p1.clone() * p3.clone();
    let y2 = mod_inverse(m2.clone() % p2.clone(), p2.clone());

    // M3 = P1 * P2, y3 = M3^(-1) mod P3
    let m3 = p1.clone() * p2.clone();
    let y3 = mod_inverse(m3.clone() % p3.clone(), p3.clone());

    // x = r1*M1*y1 + r2*M2*y2 + r3*M3*y3 (mod M)
    let term1 = Integer::from(r1 as i64) * m1 * y1;
    let term2 = Integer::from(r2 as i64) * m2 * y2;
    let term3 = Integer::from(r3 as i64) * m3 * y3;

    let x = (term1 + term2 + term3) % m.clone();

    // Adjust for signed result (assume result fits in [-M/2, M/2))
    let half_m = m.clone() / Integer::from(2i64);
    if x > half_m {
        x - m
    } else {
        x
    }
}

/// Computes the modular inverse using extended Euclidean algorithm.
fn mod_inverse(a: Integer, m: Integer) -> Integer {
    let (g, x, _) = extended_gcd(a, m.clone());
    if g != Integer::from(1i64) {
        panic!("No modular inverse exists");
    }
    // Ensure positive result
    let result = x % m.clone();
    if result.is_negative() {
        result + m
    } else {
        result
    }
}

/// Extended Euclidean algorithm: returns (gcd, x, y) such that a*x + b*y = gcd.
fn extended_gcd(a: Integer, b: Integer) -> (Integer, Integer, Integer) {
    if b.is_zero() {
        return (a, Integer::from(1i64), Integer::from(0i64));
    }

    let (g, x1, y1) = extended_gcd(b.clone(), a.clone() % b.clone());
    let x = y1.clone();
    let y = x1 - (a / b) * y1;
    (g, x, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_multiply_simple() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        let a = vec![Integer::new(1), Integer::new(2)];
        let b = vec![Integer::new(3), Integer::new(4)];

        let c = fft_multiply_integers(&a, &b);

        assert_eq!(c.len(), 3);
        assert_eq!(c[0].to_i64(), Some(3));
        assert_eq!(c[1].to_i64(), Some(10));
        assert_eq!(c[2].to_i64(), Some(8));
    }

    #[test]
    fn test_fft_multiply_negative() {
        // (1 - 2x) * (3 + 4x) = 3 - 2x - 8x^2
        let a = vec![Integer::new(1), Integer::new(-2)];
        let b = vec![Integer::new(3), Integer::new(4)];

        let c = fft_multiply_integers(&a, &b);

        assert_eq!(c.len(), 3);
        assert_eq!(c[0].to_i64(), Some(3));
        assert_eq!(c[1].to_i64(), Some(-2));
        assert_eq!(c[2].to_i64(), Some(-8));
    }

    #[test]
    fn test_fft_multiply_larger() {
        // Test with larger coefficients
        let a: Vec<Integer> = (1..=10).map(Integer::new).collect();
        let b: Vec<Integer> = (1..=10).map(Integer::new).collect();

        let c = fft_multiply_integers(&a, &b);

        // Verify against schoolbook
        let mut expected = vec![Integer::new(0); a.len() + b.len() - 1];
        for i in 0..a.len() {
            for j in 0..b.len() {
                expected[i + j] = expected[i + j].clone() + a[i].clone() * b[j].clone();
            }
        }

        assert_eq!(c.len(), expected.len());
        for (i, (got, exp)) in c.iter().zip(expected.iter()).enumerate() {
            assert_eq!(got, exp, "mismatch at index {i}");
        }
    }

    #[test]
    fn test_crt3() {
        // Test CRT reconstruction
        let x = Integer::new(12345);
        let r1 = (x.clone() % Integer::from(primes::P1 as i64))
            .to_i64()
            .unwrap() as u64;
        let r2 = (x.clone() % Integer::from(primes::P2 as i64))
            .to_i64()
            .unwrap() as u64;
        let r3 = (x.clone() % Integer::from(primes::P3 as i64))
            .to_i64()
            .unwrap() as u64;

        let reconstructed = crt3(r1, r2, r3);
        assert_eq!(reconstructed, x);
    }
}
