//! Number Theoretic Transform (NTT) for exact polynomial multiplication.
//!
//! NTT is the finite field analog of FFT, enabling exact O(n log n)
//! polynomial multiplication without floating-point errors.

use tertius_integers::ModInt;

/// NTT-friendly prime: 998244353 = 2^23 * 7 * 17 + 1
/// This prime supports NTT up to length 2^23.
pub const NTT_PRIME: u64 = 998_244_353;

/// Primitive root of NTT_PRIME.
pub const PRIMITIVE_ROOT: u64 = 3;

/// Type alias for the NTT field.
pub type NttField = ModInt<NTT_PRIME>;

/// Computes the NTT of a polynomial in-place.
///
/// The input length must be a power of 2.
/// After NTT, `a[i]` contains the evaluation at ω^i where ω is a primitive n-th root of unity.
pub fn ntt(a: &mut [NttField]) {
    let n = a.len();
    debug_assert!(n.is_power_of_two(), "NTT length must be power of 2");

    if n == 1 {
        return;
    }

    // Bit-reversal permutation
    bit_reverse(a);

    // Cooley-Tukey iterative NTT
    let mut len = 2;
    while len <= n {
        // ω_len = primitive len-th root of unity
        let w_len = NttField::new(PRIMITIVE_ROOT).pow((NTT_PRIME - 1) / len as u64);

        for i in (0..n).step_by(len) {
            let mut w = NttField::new(1);
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

/// Computes the inverse NTT in-place.
///
/// This recovers the original polynomial coefficients from the NTT representation.
pub fn intt(a: &mut [NttField]) {
    let n = a.len();
    debug_assert!(n.is_power_of_two(), "INTT length must be power of 2");

    if n == 1 {
        return;
    }

    // Bit-reversal permutation
    bit_reverse(a);

    // Cooley-Tukey iterative INTT (using inverse root)
    let mut len = 2;
    while len <= n {
        // ω_len^(-1) = inverse of primitive len-th root of unity
        let w_len = NttField::new(PRIMITIVE_ROOT)
            .pow((NTT_PRIME - 1) / len as u64)
            .inv()
            .unwrap();

        for i in (0..n).step_by(len) {
            let mut w = NttField::new(1);
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
    let n_inv = NttField::new(n as u64).inv().unwrap();
    for x in a.iter_mut() {
        *x = *x * n_inv;
    }
}

/// Performs bit-reversal permutation in-place.
fn bit_reverse(a: &mut [NttField]) {
    let n = a.len();
    let log_n = n.trailing_zeros();

    for i in 0..n {
        let j = reverse_bits(i as u32, log_n) as usize;
        if i < j {
            a.swap(i, j);
        }
    }
}

/// Reverses the lower `bits` bits of `x`.
#[inline]
fn reverse_bits(x: u32, bits: u32) -> u32 {
    x.reverse_bits() >> (32 - bits)
}

/// Multiplies two polynomials using NTT.
///
/// Returns coefficients modulo NTT_PRIME.
pub fn ntt_multiply(a: &[NttField], b: &[NttField]) -> Vec<NttField> {
    if a.is_empty() || b.is_empty() {
        return vec![NttField::new(0)];
    }

    // Result length
    let result_len = a.len() + b.len() - 1;
    let n = result_len.next_power_of_two();

    // Pad to power of 2
    let mut a_padded = Vec::with_capacity(n);
    a_padded.extend_from_slice(a);
    a_padded.resize(n, NttField::new(0));

    let mut b_padded = Vec::with_capacity(n);
    b_padded.extend_from_slice(b);
    b_padded.resize(n, NttField::new(0));

    // Forward NTT
    ntt(&mut a_padded);
    ntt(&mut b_padded);

    // Pointwise multiplication
    for i in 0..n {
        a_padded[i] = a_padded[i] * b_padded[i];
    }

    // Inverse NTT
    intt(&mut a_padded);

    // Trim to actual result length
    a_padded.truncate(result_len);
    a_padded
}

/// Additional NTT-friendly primes for CRT-based multiplication.
pub mod primes {
    /// 2^23 * 119 + 1 = 998244353 (primary)
    pub const P1: u64 = 998_244_353;
    /// 2^24 * 73 + 1 = 1224736769
    pub const P2: u64 = 1_224_736_769;
    /// 2^26 * 7 + 1 = 469762049
    pub const P3: u64 = 469_762_049;

    /// Primitive root for P1.
    pub const G1: u64 = 3;
    /// Primitive root for P2.
    pub const G2: u64 = 3;
    /// Primitive root for P3.
    pub const G3: u64 = 3;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_intt_roundtrip() {
        let original: Vec<NttField> = (0..8)
            .map(|i| NttField::new(i + 1))
            .collect();

        let mut a = original.clone();
        ntt(&mut a);
        intt(&mut a);

        for (i, (orig, recovered)) in original.iter().zip(a.iter()).enumerate() {
            assert_eq!(orig.value(), recovered.value(), "mismatch at index {i}");
        }
    }

    #[test]
    fn test_ntt_multiply_simple() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        let a = vec![NttField::new(1), NttField::new(2)];
        let b = vec![NttField::new(3), NttField::new(4)];

        let c = ntt_multiply(&a, &b);

        assert_eq!(c.len(), 3);
        assert_eq!(c[0].value(), 3);
        assert_eq!(c[1].value(), 10);
        assert_eq!(c[2].value(), 8);
    }

    #[test]
    fn test_ntt_multiply_larger() {
        // (1 + x + x^2 + x^3) * (1 + x + x^2 + x^3)
        // = 1 + 2x + 3x^2 + 4x^3 + 3x^4 + 2x^5 + x^6
        let a: Vec<NttField> = (0..4).map(|_| NttField::new(1)).collect();
        let b = a.clone();

        let c = ntt_multiply(&a, &b);

        assert_eq!(c.len(), 7);
        assert_eq!(c[0].value(), 1);
        assert_eq!(c[1].value(), 2);
        assert_eq!(c[2].value(), 3);
        assert_eq!(c[3].value(), 4);
        assert_eq!(c[4].value(), 3);
        assert_eq!(c[5].value(), 2);
        assert_eq!(c[6].value(), 1);
    }

    #[test]
    fn test_bit_reverse() {
        let mut a: Vec<NttField> = (0..8).map(|i| NttField::new(i)).collect();
        bit_reverse(&mut a);

        // For n=8, bit reversal is: 0->0, 1->4, 2->2, 3->6, 4->1, 5->5, 6->3, 7->7
        let expected = [0, 4, 2, 6, 1, 5, 3, 7];
        for (i, &exp) in expected.iter().enumerate() {
            assert_eq!(a[i].value(), exp, "mismatch at index {i}");
        }
    }
}
