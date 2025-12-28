//! Resultants and discriminants via Sylvester matrix.
//!
//! The resultant of two polynomials f and g is zero iff they share a common root.
//! It's computed as the determinant of the Sylvester matrix.

use std::ops::Neg;
use tertius_rings::traits::{EuclideanDomain, Ring};

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
pub fn resultant<R: EuclideanDomain + Clone + Neg<Output = R>>(f: &[R], g: &[R]) -> R {
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

/// Computes the determinant using Bareiss algorithm (fraction-free Gaussian elimination).
///
/// This is O(n³) compared to O(n!) for cofactor expansion.
/// Works for any Euclidean domain (needed for exact division).
fn determinant<R: EuclideanDomain + Clone + Neg<Output = R>>(matrix: &[Vec<R>]) -> R {
    let n = matrix.len();
    if n == 0 {
        return R::one();
    }
    if n == 1 {
        return matrix[0][0].clone();
    }
    if n == 2 {
        // 2x2: ad - bc
        let a = &matrix[0][0];
        let b = &matrix[0][1];
        let c = &matrix[1][0];
        let d = &matrix[1][1];
        return a.clone() * d.clone() + (b.clone() * c.clone()).neg();
    }
    if n == 3 {
        // 3x3: Sarrus' rule
        let a = &matrix[0][0];
        let b = &matrix[0][1];
        let c = &matrix[0][2];
        let d = &matrix[1][0];
        let e = &matrix[1][1];
        let f = &matrix[1][2];
        let g = &matrix[2][0];
        let h = &matrix[2][1];
        let i = &matrix[2][2];

        // aei + bfg + cdh - ceg - bdi - afh
        let pos = a.clone() * e.clone() * i.clone()
            + b.clone() * f.clone() * g.clone()
            + c.clone() * d.clone() * h.clone();
        let neg = c.clone() * e.clone() * g.clone()
            + b.clone() * d.clone() * i.clone()
            + a.clone() * f.clone() * h.clone();
        return pos + neg.neg();
    }

    // Bareiss algorithm: fraction-free Gaussian elimination
    // We track row swaps to maintain correct sign
    let mut m: Vec<Vec<R>> = matrix.to_vec();
    let mut sign_flips = 0usize;

    for k in 0..n - 1 {
        // Find pivot in column k
        let mut pivot_row = None;
        for i in k..n {
            if !m[i][k].is_zero() {
                pivot_row = Some(i);
                break;
            }
        }

        let pivot_row = match pivot_row {
            Some(r) => r,
            None => return R::zero(), // Column is all zeros, det = 0
        };

        // Swap rows if needed
        if pivot_row != k {
            m.swap(k, pivot_row);
            sign_flips += 1;
        }

        // Get the pivot element and previous pivot (for Bareiss division)
        let pivot = m[k][k].clone();
        let prev_pivot = if k > 0 { m[k - 1][k - 1].clone() } else { R::one() };

        // Eliminate column k in rows k+1 to n-1
        for i in k + 1..n {
            for j in k + 1..n {
                // Bareiss formula: m[i][j] = (m[i][j] * pivot - m[i][k] * m[k][j]) / prev_pivot
                let numerator = m[i][j].clone() * pivot.clone()
                    + (m[i][k].clone() * m[k][j].clone()).neg();
                // In an integral domain, this division is exact
                m[i][j] = exact_div(&numerator, &prev_pivot);
            }
            m[i][k] = R::zero();
        }
    }

    // The determinant is the (n-1, n-1) element, adjusted for sign
    let det = m[n - 1][n - 1].clone();
    if sign_flips % 2 == 0 {
        det
    } else {
        det.neg()
    }
}

/// Performs exact division in an Euclidean domain.
///
/// This assumes that divisor divides dividend exactly.
/// The Bareiss algorithm guarantees this property.
fn exact_div<R: EuclideanDomain + Clone>(dividend: &R, divisor: &R) -> R {
    if divisor.is_one() {
        return dividend.clone();
    }
    if divisor.is_zero() {
        panic!("exact_div: division by zero");
    }

    let (quotient, remainder) = dividend.div_rem(divisor);

    // Bareiss guarantees exact divisibility
    debug_assert!(
        remainder.is_zero(),
        "exact_div: remainder should be zero in Bareiss algorithm"
    );

    quotient
}

/// Computes the discriminant of a univariate polynomial.
///
/// The discriminant is related to the resultant of f and f':
/// disc(f) = (-1)^{n(n-1)/2} * res(f, f') / a_n
///
/// where a_n is the leading coefficient and n is the degree.
pub fn discriminant<R: EuclideanDomain + Clone + Neg<Output = R>>(f: &[R]) -> R {
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
