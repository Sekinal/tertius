//! Berlekamp-Massey algorithm for finding minimal polynomials.
//!
//! Given a sequence s_0, s_1, ..., s_{2n-1}, this algorithm finds the
//! minimal polynomial C(x) = 1 + c_1*x + ... + c_L*x^L such that
//! s_i = -sum_{j=1}^L c_j * s_{i-j} for all i >= L.
//!
//! Used in:
//! - Block Wiedemann algorithm
//! - Ben-Or/Tiwari sparse interpolation
//! - Decoding BCH and Reed-Solomon codes

use std::ops::Neg;
use tertius_rings::traits::{Field, Ring};

/// Result of the Berlekamp-Massey algorithm.
#[derive(Clone, Debug)]
pub struct BerlekampMasseyResult<R> {
    /// The connection polynomial C(x) = c_0 + c_1*x + ... + c_L*x^L.
    /// c_0 = 1 always.
    pub connection_poly: Vec<R>,
    /// The linear complexity L.
    pub linear_complexity: usize,
}

impl<R: Clone> BerlekampMasseyResult<R> {
    /// Returns the coefficients of the connection polynomial.
    pub fn coeffs(&self) -> &[R] {
        &self.connection_poly
    }

    /// Returns the degree of the connection polynomial.
    pub fn degree(&self) -> usize {
        self.linear_complexity
    }
}

/// Runs the Berlekamp-Massey algorithm on a sequence over a field.
///
/// # Arguments
/// - `sequence`: The sequence s_0, s_1, ..., s_{2n-1}
///
/// # Returns
/// The minimal connection polynomial and its degree (linear complexity).
pub fn berlekamp_massey<R>(sequence: &[R]) -> BerlekampMasseyResult<R>
where
    R: Field + Clone + Neg<Output = R>,
{
    let n = sequence.len();

    if n == 0 {
        return BerlekampMasseyResult {
            connection_poly: vec![R::one()],
            linear_complexity: 0,
        };
    }

    // C(x) = current connection polynomial
    let mut c: Vec<R> = vec![R::one()];
    // B(x) = previous connection polynomial before last length change
    let mut b: Vec<R> = vec![R::one()];
    // L = current linear complexity
    let mut l = 0usize;
    // m = number of iterations since last length change
    let mut m = 1usize;
    // d = last non-zero discrepancy
    let mut d_b = R::one();

    for i in 0..n {
        // Compute discrepancy d = s_i + sum_{j=1}^L c_j * s_{i-j}
        let mut d = sequence[i].clone();
        for j in 1..c.len().min(i + 1) {
            d = d + c[j].clone() * sequence[i - j].clone();
        }

        if d.is_zero() {
            // No update needed
            m += 1;
        } else if 2 * l <= i {
            // Need to increase complexity
            let t = c.clone();

            // C(x) = C(x) - (d / d_b) * x^m * B(x)
            let scale = d.clone() * d_b.inv().expect("d_b must be invertible");
            let neg_scale = scale.neg();

            // Ensure c is long enough
            while c.len() < b.len() + m {
                c.push(R::zero());
            }

            // Subtract scale * x^m * B(x) from C(x)
            for (j, bj) in b.iter().enumerate() {
                c[j + m] = c[j + m].clone() + neg_scale.clone() * bj.clone();
            }

            l = i + 1 - l;
            b = t;
            d_b = d;
            m = 1;
        } else {
            // C(x) = C(x) - (d / d_b) * x^m * B(x)
            let scale = d.clone() * d_b.inv().expect("d_b must be invertible");
            let neg_scale = scale.neg();

            // Ensure c is long enough
            while c.len() < b.len() + m {
                c.push(R::zero());
            }

            // Subtract scale * x^m * B(x) from C(x)
            for (j, bj) in b.iter().enumerate() {
                c[j + m] = c[j + m].clone() + neg_scale.clone() * bj.clone();
            }

            m += 1;
        }
    }

    // Trim trailing zeros
    while c.len() > 1 && c.last().map_or(false, |x| x.is_zero()) {
        c.pop();
    }

    BerlekampMasseyResult {
        connection_poly: c,
        linear_complexity: l,
    }
}

/// Finds the roots of a polynomial over a finite field using exhaustive search.
///
/// This is simple but O(n * |F|). For large fields, use Cantor-Zassenhaus.
pub fn find_roots_exhaustive<const P: u64>(
    poly: &[tertius_rings::finite_field::FiniteField<P>],
) -> Vec<tertius_rings::finite_field::FiniteField<P>> {
    use tertius_rings::finite_field::FiniteField;

    let mut roots = Vec::new();

    for val in 0..P {
        let x = FiniteField::new(val);

        // Evaluate polynomial at x using Horner's method
        let mut result = FiniteField::<P>::zero();
        for coeff in poly.iter().rev() {
            result = result * x.clone() + coeff.clone();
        }

        if result.is_zero() {
            roots.push(x);
        }
    }

    roots
}

/// Reconstructs the sequence elements from a connection polynomial.
///
/// Given C(x) and initial values, computes more elements of the sequence.
pub fn extend_sequence<R>(
    connection_poly: &[R],
    initial: &[R],
    target_length: usize,
) -> Vec<R>
where
    R: Ring + Clone + Neg<Output = R>,
{
    let l = connection_poly.len() - 1; // Degree of C(x)

    if initial.len() < l {
        panic!("Need at least {} initial values", l);
    }

    let mut sequence = initial.to_vec();

    while sequence.len() < target_length {
        let i = sequence.len();

        // s_i = -sum_{j=1}^L c_j * s_{i-j}
        let mut val = R::zero();
        for j in 1..=l {
            val = val + connection_poly[j].clone() * sequence[i - j].clone();
        }
        sequence.push(val.neg());
    }

    sequence
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::finite_field::FiniteField;

    type GF7 = FiniteField<7>;

    fn gf7(n: u64) -> GF7 {
        GF7::new(n % 7)
    }

    #[test]
    fn test_berlekamp_massey_fibonacci() {
        // Fibonacci over GF(7): 1, 1, 2, 3, 5, 1, 6, 0, 6, 6, ...
        // Connection polynomial: C(x) = 1 - x - x^2
        let fib: Vec<GF7> = vec![
            gf7(1),
            gf7(1),
            gf7(2),
            gf7(3),
            gf7(5),
            gf7(1), // 8 mod 7
            gf7(6), // 13 mod 7
            gf7(0), // 21 mod 7
        ];

        let result = berlekamp_massey(&fib);

        // C(x) = 1 - x - x^2
        assert_eq!(result.linear_complexity, 2);
        assert_eq!(result.connection_poly.len(), 3);
        assert_eq!(result.connection_poly[0], gf7(1));
        assert_eq!(result.connection_poly[1], gf7(6)); // -1 mod 7
        assert_eq!(result.connection_poly[2], gf7(6)); // -1 mod 7
    }

    #[test]
    fn test_berlekamp_massey_constant() {
        // Constant sequence: 3, 3, 3, 3
        // Connection polynomial: C(x) = 1 - x
        let seq: Vec<GF7> = vec![gf7(3), gf7(3), gf7(3), gf7(3)];

        let result = berlekamp_massey(&seq);

        assert_eq!(result.linear_complexity, 1);
    }

    #[test]
    fn test_berlekamp_massey_alternating() {
        // Alternating: 1, -1, 1, -1 = 1, 6, 1, 6 in GF(7)
        // Connection polynomial: C(x) = 1 + x
        let seq: Vec<GF7> = vec![gf7(1), gf7(6), gf7(1), gf7(6)];

        let result = berlekamp_massey(&seq);

        assert_eq!(result.linear_complexity, 1);
    }

    #[test]
    fn test_extend_sequence() {
        // Fibonacci: s_n = s_{n-1} + s_{n-2}
        // C(x) = 1 - x - x^2
        let conn: Vec<GF7> = vec![gf7(1), gf7(6), gf7(6)];
        let initial: Vec<GF7> = vec![gf7(1), gf7(1)];

        let extended = extend_sequence(&conn, &initial, 8);

        assert_eq!(extended[0], gf7(1));
        assert_eq!(extended[1], gf7(1));
        assert_eq!(extended[2], gf7(2));
        assert_eq!(extended[3], gf7(3));
        assert_eq!(extended[4], gf7(5));
    }
}
