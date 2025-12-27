//! Karatsuba multiplication algorithm.
//!
//! This module provides the Karatsuba divide-and-conquer multiplication
//! algorithm, which achieves O(n^1.58) complexity.

use tertius_rings::traits::Ring;

/// Karatsuba multiplication threshold.
///
/// Below this degree, schoolbook multiplication is faster.
pub const KARATSUBA_THRESHOLD: usize = 32;

/// Performs Karatsuba multiplication on coefficient slices.
///
/// Returns a new vector containing the product coefficients.
pub fn karatsuba_mul<R: Ring>(a: &[R], b: &[R]) -> Vec<R> {
    let n = a.len();
    let m = b.len();

    // Base case
    if n < KARATSUBA_THRESHOLD || m < KARATSUBA_THRESHOLD {
        return schoolbook_mul(a, b);
    }

    // Make both the same size (power of 2)
    let size = n.max(m).next_power_of_two();
    let half = size / 2;

    // Extend a and b to size
    let mut a_ext = a.to_vec();
    let mut b_ext = b.to_vec();
    a_ext.resize(size, R::zero());
    b_ext.resize(size, R::zero());

    // Split: a = a0 + a1*x^half, b = b0 + b1*x^half
    let a0 = &a_ext[..half];
    let a1 = &a_ext[half..];
    let b0 = &b_ext[..half];
    let b1 = &b_ext[half..];

    // Compute z0 = a0*b0, z2 = a1*b1
    let z0 = karatsuba_mul(a0, b0);
    let z2 = karatsuba_mul(a1, b1);

    // Compute (a0+a1) and (b0+b1)
    let a01: Vec<R> = a0.iter().zip(a1.iter())
        .map(|(x, y)| x.clone() + y.clone())
        .collect();
    let b01: Vec<R> = b0.iter().zip(b1.iter())
        .map(|(x, y)| x.clone() + y.clone())
        .collect();

    // z1 = (a0+a1)*(b0+b1) - z0 - z2
    let z1_sum = karatsuba_mul(&a01, &b01);
    let mut z1 = z1_sum;

    // Subtract z0 and z2 from z1
    for (i, c) in z0.iter().enumerate() {
        if i < z1.len() {
            z1[i] = z1[i].clone() - c.clone();
        }
    }
    for (i, c) in z2.iter().enumerate() {
        if i < z1.len() {
            z1[i] = z1[i].clone() - c.clone();
        }
    }

    // Combine: result = z0 + z1*x^half + z2*x^(2*half)
    let result_len = 2 * size - 1;
    let mut result = vec![R::zero(); result_len];

    for (i, c) in z0.into_iter().enumerate() {
        result[i] = c;
    }

    for (i, c) in z1.into_iter().enumerate() {
        result[i + half] = result[i + half].clone() + c;
    }

    for (i, c) in z2.into_iter().enumerate() {
        result[i + 2 * half] = result[i + 2 * half].clone() + c;
    }

    // Trim trailing zeros
    while result.len() > 1 && result.last().map_or(false, |c| c.is_zero()) {
        result.pop();
    }

    result
}

/// Schoolbook multiplication: O(n²).
pub fn schoolbook_mul<R: Ring>(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return vec![R::zero()];
    }

    let n = a.len();
    let m = b.len();
    let mut result = vec![R::zero(); n + m - 1];

    for i in 0..n {
        for j in 0..m {
            result[i + j] = result[i + j].clone() + a[i].clone() * b[j].clone();
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_schoolbook() {
        let a = vec![Q::from_integer(1), Q::from_integer(2)]; // 1 + 2x
        let b = vec![Q::from_integer(3), Q::from_integer(4)]; // 3 + 4x

        let c = schoolbook_mul(&a, &b);
        // (1 + 2x)(3 + 4x) = 3 + 10x + 8x^2
        assert_eq!(c[0], Q::from_integer(3));
        assert_eq!(c[1], Q::from_integer(10));
        assert_eq!(c[2], Q::from_integer(8));
    }

    #[test]
    fn test_karatsuba_small() {
        // Small polynomials should still work
        let a = vec![Q::from_integer(1), Q::from_integer(2)];
        let b = vec![Q::from_integer(3), Q::from_integer(4)];

        let c = karatsuba_mul(&a, &b);
        assert_eq!(c[0], Q::from_integer(3));
        assert_eq!(c[1], Q::from_integer(10));
        assert_eq!(c[2], Q::from_integer(8));
    }

    #[test]
    fn test_karatsuba_large() {
        // Create larger polynomials to trigger Karatsuba
        let n = 100;
        let a: Vec<Q> = (0..n).map(|i| Q::from_integer(i)).collect();
        let b: Vec<Q> = (0..n).map(|i| Q::from_integer(n - i)).collect();

        let c_school = schoolbook_mul(&a, &b);
        let c_kara = karatsuba_mul(&a, &b);

        // Results should match
        assert_eq!(c_school.len(), c_kara.len());
        for i in 0..c_school.len() {
            assert_eq!(c_school[i], c_kara[i], "mismatch at index {i}");
        }
    }
}
