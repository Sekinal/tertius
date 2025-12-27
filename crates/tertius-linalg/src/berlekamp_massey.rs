//! Berlekamp-Massey algorithm for finding minimal polynomials.
//!
//! Given a sequence s_0, s_1, ..., s_{n-1}, this algorithm finds the
//! minimal polynomial C(x) such that:
//!   C_0 * s_j + C_1 * s_{j-1} + ... + C_L * s_{j-L} = 0
//! for all j >= L.
//!
//! This is used by:
//! - Block Wiedemann algorithm for sparse linear algebra
//! - Ben-Or/Tiwari sparse polynomial interpolation

use tertius_rings::traits::Field;

/// Result of the Berlekamp-Massey algorithm.
#[derive(Clone, Debug)]
pub struct BerlekampMasseyResult<R> {
    /// Connection polynomial C(x) = C_0 + C_1*x + ... + C_L*x^L.
    /// Coefficients in ascending degree order.
    pub connection_poly: Vec<R>,
    /// Length of the linear recurrence (degree of connection polynomial).
    pub length: usize,
}

/// Computes the minimal polynomial of a linear recurrence sequence.
///
/// Given a sequence `s = [s_0, s_1, ..., s_{n-1}]`, finds the shortest
/// linear recurrence relation satisfied by the sequence.
///
/// # Returns
///
/// The connection polynomial C(x) such that:
///   C_0 * s_j + C_1 * s_{j-1} + ... + C_L * s_{j-L} = 0
///
/// The polynomial is returned with coefficients in ascending degree order.
#[must_use]
pub fn berlekamp_massey<R: Field + Clone>(sequence: &[R]) -> BerlekampMasseyResult<R> {
    if sequence.is_empty() {
        return BerlekampMasseyResult {
            connection_poly: vec![R::one()],
            length: 0,
        };
    }

    let n = sequence.len();

    // Current connection polynomial C(x)
    let mut c = vec![R::one()];

    // Previous connection polynomial B(x)
    let mut b = vec![R::one()];

    // Current recurrence length
    let mut l = 0usize;

    // Number of iterations since last length change
    let mut m = 1usize;

    // Previous discrepancy
    let mut delta_prev = R::one();

    for i in 0..n {
        // Compute discrepancy
        let delta = compute_discrepancy(&c, sequence, i);

        if delta.is_zero() {
            m += 1;
        } else if 2 * l <= i {
            // Length needs to increase
            let t = c.clone();

            // C(x) = C(x) - (delta / delta_prev) * x^m * B(x)
            let scale = delta.field_div(&delta_prev);
            c = poly_sub_scaled_shift(&c, &b, &scale, m);

            // Update B and other state
            b = t;
            l = i + 1 - l;
            delta_prev = delta;
            m = 1;
        } else {
            // Just update C(x)
            let scale = delta.field_div(&delta_prev);
            c = poly_sub_scaled_shift(&c, &b, &scale, m);
            m += 1;
        }
    }

    // Normalize so that leading coefficient is 1 (monic)
    normalize_monic(&mut c);

    BerlekampMasseyResult {
        length: l,
        connection_poly: c,
    }
}

/// Computes the discrepancy at step n.
fn compute_discrepancy<R: Field + Clone>(c: &[R], s: &[R], n: usize) -> R {
    let mut delta = R::zero();

    for (j, cj) in c.iter().enumerate() {
        if j <= n {
            delta = delta + cj.clone() * s[n - j].clone();
        }
    }

    delta
}

/// Computes A(x) - scale * x^shift * B(x).
fn poly_sub_scaled_shift<R: Field + Clone>(a: &[R], b: &[R], scale: &R, shift: usize) -> Vec<R> {
    let result_len = a.len().max(b.len() + shift);
    let mut result = vec![R::zero(); result_len];

    // Copy A
    for (i, ai) in a.iter().enumerate() {
        result[i] = ai.clone();
    }

    // Subtract scale * x^shift * B
    for (i, bi) in b.iter().enumerate() {
        let idx = i + shift;
        if idx < result.len() {
            result[idx] = result[idx].clone() + (bi.clone() * scale.clone()).neg();
        }
    }

    // Trim trailing zeros
    while result.len() > 1 && result.last().map_or(false, |x| x.is_zero()) {
        result.pop();
    }

    result
}

/// Normalizes a polynomial to be monic (leading coefficient = 1).
fn normalize_monic<R: Field + Clone>(poly: &mut Vec<R>) {
    if poly.is_empty() {
        return;
    }

    // Find last non-zero coefficient
    while poly.len() > 1 && poly.last().map_or(false, |x| x.is_zero()) {
        poly.pop();
    }

    if let Some(lead) = poly.last().cloned() {
        if !lead.is_zero() {
            if let Some(inv) = lead.inv() {
                for coeff in poly.iter_mut() {
                    *coeff = coeff.clone() * inv.clone();
                }
            }
        }
    }
}

/// Matrix version of Berlekamp-Massey for block algorithms.
///
/// Given a sequence of matrices S_0, S_1, ..., finds the minimal
/// matrix polynomial.
#[derive(Clone, Debug)]
pub struct MatrixBMResult<R> {
    /// Matrix connection polynomial coefficients.
    pub coefficients: Vec<Vec<Vec<R>>>,
    /// Degree of the polynomial.
    pub degree: usize,
}

/// Berlekamp-Massey for matrix sequences (simplified version).
///
/// Used in Block Wiedemann when we have block_size x block_size matrices.
#[must_use]
pub fn matrix_berlekamp_massey<R: Field + Clone>(
    sequence: &[Vec<Vec<R>>],
    block_size: usize,
) -> MatrixBMResult<R> {
    // For the full matrix BM algorithm, we would need to track
    // separate polynomials for each row. This is a simplified version
    // that processes each row independently.

    if sequence.is_empty() || block_size == 0 {
        return MatrixBMResult {
            coefficients: vec![],
            degree: 0,
        };
    }

    // Extract the first row of each matrix and apply scalar BM
    // This is a simplification - full matrix BM is more complex
    let first_row_seq: Vec<R> = sequence
        .iter()
        .filter_map(|mat| mat.first().and_then(|row| row.first().cloned()))
        .collect();

    let result = berlekamp_massey(&first_row_seq);

    // Convert to matrix form (diagonal matrices with scalar polynomial)
    let coefficients: Vec<Vec<Vec<R>>> = result
        .connection_poly
        .iter()
        .map(|c| {
            let mut mat = vec![vec![R::zero(); block_size]; block_size];
            for i in 0..block_size {
                mat[i][i] = c.clone();
            }
            mat
        })
        .collect();

    MatrixBMResult {
        degree: result.length,
        coefficients,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;
    use tertius_rings::traits::Ring;

    #[test]
    fn test_constant_sequence() {
        // Sequence [1, 1, 1, 1] has minimal polynomial C(x) = 1 - x
        // (with convention C_0 * s_n + C_1 * s_{n-1} = 0)
        let seq: Vec<Q> = vec![
            Q::from_integer(1),
            Q::from_integer(1),
            Q::from_integer(1),
            Q::from_integer(1),
        ];

        let result = berlekamp_massey(&seq);
        assert_eq!(result.length, 1);
    }

    #[test]
    fn test_fibonacci_like() {
        // Fibonacci: s_n = s_{n-1} + s_{n-2}
        // Connection polynomial: 1 - x - x^2
        let seq: Vec<Q> = vec![
            Q::from_integer(1),
            Q::from_integer(1),
            Q::from_integer(2),
            Q::from_integer(3),
            Q::from_integer(5),
            Q::from_integer(8),
        ];

        let result = berlekamp_massey(&seq);
        assert_eq!(result.length, 2);
    }

    #[test]
    fn test_geometric() {
        // Geometric sequence s_n = 2^n: s_n = 2 * s_{n-1}
        // Connection polynomial: 1 - 2x
        let seq: Vec<Q> = vec![
            Q::from_integer(1),
            Q::from_integer(2),
            Q::from_integer(4),
            Q::from_integer(8),
            Q::from_integer(16),
        ];

        let result = berlekamp_massey(&seq);
        assert_eq!(result.length, 1);

        // Verify connection polynomial
        assert_eq!(result.connection_poly.len(), 2);
    }

    #[test]
    fn test_empty_sequence() {
        let seq: Vec<Q> = vec![];
        let result = berlekamp_massey(&seq);
        assert_eq!(result.length, 0);
        assert_eq!(result.connection_poly.len(), 1);
    }

    #[test]
    fn test_verify_recurrence() {
        // Verify that the computed polynomial satisfies the recurrence
        let seq: Vec<Q> = vec![
            Q::from_integer(1),
            Q::from_integer(3),
            Q::from_integer(7),
            Q::from_integer(15),
            Q::from_integer(31),
            Q::from_integer(63),
        ];

        let result = berlekamp_massey(&seq);

        // Check that C_0 * s_j + C_1 * s_{j-1} + ... = 0 for j >= L
        for j in result.length..seq.len() {
            let mut sum = Q::from_integer(0);
            for (k, ck) in result.connection_poly.iter().enumerate() {
                if k <= j {
                    sum = sum + ck.clone() * seq[j - k].clone();
                }
            }
            assert!(sum.is_zero(), "Recurrence not satisfied at j={}", j);
        }
    }
}
