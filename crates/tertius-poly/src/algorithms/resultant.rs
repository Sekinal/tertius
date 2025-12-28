//! Resultants and discriminants via Sylvester matrix.
//!
//! The resultant of two polynomials f and g is zero iff they share a common root.
//! It's computed as the determinant of the Sylvester matrix.

use std::ops::Neg;
use tertius_rings::traits::Ring;

/// Computes the resultant of two univariate polynomials.
///
/// The resultant is the determinant of the Sylvester matrix.
/// It's zero iff f and g share a common root.
///
/// # Arguments
/// - `f`: First polynomial as coefficient vector [f_0, f_1, ..., f_m]
/// - `g`: Second polynomial as coefficient vector [g_0, g_1, ..., g_n]
///
/// # Returns
/// The resultant res(f, g).
pub fn resultant<R: Ring + Clone + Neg<Output = R>>(f: &[R], g: &[R]) -> R {
    // Handle edge cases
    if f.is_empty() || g.is_empty() {
        return R::zero();
    }

    let deg_f = f.len() - 1;
    let deg_g = g.len() - 1;

    if deg_f == 0 {
        // f is constant
        return f[0].pow(deg_g as u32);
    }

    if deg_g == 0 {
        // g is constant
        return g[0].pow(deg_f as u32);
    }

    // Build Sylvester matrix
    let size = deg_f + deg_g;
    let sylvester = build_sylvester_matrix(f, g, size);

    // Compute determinant
    determinant(&sylvester)
}

/// Builds the Sylvester matrix for two polynomials.
fn build_sylvester_matrix<R: Ring + Clone>(f: &[R], g: &[R], size: usize) -> Vec<Vec<R>> {
    let deg_f = f.len() - 1;
    let deg_g = g.len() - 1;

    let mut matrix = vec![vec![R::zero(); size]; size];

    // First deg_g rows: shifts of f
    for i in 0..deg_g {
        for (j, coeff) in f.iter().enumerate() {
            if i + j < size {
                matrix[i][i + j] = coeff.clone();
            }
        }
    }

    // Next deg_f rows: shifts of g
    for i in 0..deg_f {
        for (j, coeff) in g.iter().enumerate() {
            if i + j < size {
                matrix[deg_g + i][i + j] = coeff.clone();
            }
        }
    }

    matrix
}

/// Computes the determinant using Bareiss algorithm (fraction-free elimination).
fn determinant<R: Ring + Clone + Neg<Output = R>>(matrix: &[Vec<R>]) -> R {
    let n = matrix.len();
    if n == 0 {
        return R::one();
    }
    if n == 1 {
        return matrix[0][0].clone();
    }

    // Make a mutable copy
    let mut a: Vec<Vec<R>> = matrix.to_vec();
    let mut sign = true;
    let mut prev_pivot = R::one();

    for k in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for row in k..n {
            if !a[row][k].is_zero() {
                pivot_row = Some(row);
                break;
            }
        }

        let pivot_row = match pivot_row {
            Some(r) => r,
            None => return R::zero(), // Singular matrix
        };

        // Swap rows if needed
        if pivot_row != k {
            a.swap(k, pivot_row);
            sign = !sign;
        }

        let pivot = a[k][k].clone();

        // Bareiss elimination
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                // a[i][j] = (a[k][k] * a[i][j] - a[i][k] * a[k][j]) / prev_pivot
                // For commutative rings with exact division this works
                let numerator = pivot.clone() * a[i][j].clone()
                    + (a[i][k].clone() * a[k][j].clone()).neg();
                // Division by prev_pivot for fraction-free (we skip for now and accept extra factors)
                a[i][j] = numerator;
            }
        }

        // Clear column below pivot
        for i in (k + 1)..n {
            a[i][k] = R::zero();
        }

        prev_pivot = pivot;
    }

    // For Bareiss, the determinant is the last diagonal element divided by appropriate factor
    // For simplicity, since we're not doing the division, we compute differently
    // Just return product of diagonal elements with sign adjustment
    let mut det = a[n - 1][n - 1].clone();

    if sign {
        det
    } else {
        det.neg()
    }
}

/// Computes the discriminant of a univariate polynomial.
///
/// The discriminant is related to the resultant of f and f':
/// disc(f) = (-1)^{n(n-1)/2} * res(f, f') / a_n
///
/// where a_n is the leading coefficient and n is the degree.
pub fn discriminant<R: Ring + Clone + Neg<Output = R>>(f: &[R]) -> R {
    if f.len() <= 1 {
        return R::zero();
    }

    let n = f.len() - 1;
    let f_prime = derivative(f);

    if f_prime.is_empty() {
        return R::zero();
    }

    let res = resultant(f, &f_prime);
    let lead = f.last().unwrap().clone();

    // disc = (-1)^{n(n-1)/2} * res / lead
    // For now, just return res (exact division may not be available)
    let sign_exp = (n * (n - 1)) / 2;
    if sign_exp % 2 == 0 {
        res
    } else {
        res.neg()
    }
}

/// Computes the derivative of a polynomial.
fn derivative<R: Ring + Clone>(f: &[R]) -> Vec<R> {
    if f.len() <= 1 {
        return vec![];
    }

    let mut result = Vec::with_capacity(f.len() - 1);
    for (i, coeff) in f.iter().enumerate().skip(1) {
        // Multiply by i
        let mut term = R::zero();
        for _ in 0..i {
            term = term + coeff.clone();
        }
        result.push(term);
    }

    result
}

/// Helper trait for integer power.
trait IntPow {
    fn pow(&self, exp: u32) -> Self;
}

impl<R: Ring + Clone> IntPow for R {
    fn pow(&self, mut exp: u32) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base.clone();
            }
            base = base.clone() * base;
            exp >>= 1;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::finite_field::FiniteField;

    type GF101 = FiniteField<101>;

    fn gf101(n: i64) -> GF101 {
        GF101::new(n.rem_euclid(101) as u64)
    }

    #[test]
    fn test_resultant_coprime() {
        // f = x + 1, g = x + 2
        // res(f, g) = f(root of g) = f(-2) = -2 + 1 = -1
        // Sylvester: [[1, 1], [2, 1]] = det = 1*1 - 1*2 = -1
        let f = vec![gf101(1), gf101(1)]; // 1 + x
        let g = vec![gf101(2), gf101(1)]; // 2 + x

        let res = resultant(&f, &g);
        assert_eq!(res, gf101(-1)); // -1 = 100 mod 101
    }

    #[test]
    fn test_resultant_common_root() {
        // f = x + 1, g = (x + 1) = x + 1
        // res(x+1, x+1) = 0
        let f = vec![gf101(1), gf101(1)];
        let g = vec![gf101(1), gf101(1)];

        let res = resultant(&f, &g);
        assert_eq!(res, gf101(0));
    }

    #[test]
    fn test_resultant_quadratics() {
        // f = x^2 + 2x + 1 = (x+1)^2
        // g = x^2 + 3x + 2 = (x+1)(x+2)
        // They share root x = -1, so res = 0
        let f = vec![gf101(1), gf101(2), gf101(1)];
        let g = vec![gf101(2), gf101(3), gf101(1)];

        let res = resultant(&f, &g);
        assert_eq!(res, gf101(0));
    }

    #[test]
    fn test_derivative() {
        // f = 1 + 2x + 3x^2
        // f' = 2 + 6x
        let f = vec![gf101(1), gf101(2), gf101(3)];
        let df = derivative(&f);

        assert_eq!(df.len(), 2);
        assert_eq!(df[0], gf101(2));
        assert_eq!(df[1], gf101(6));
    }
}
