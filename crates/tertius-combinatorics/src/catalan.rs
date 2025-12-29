//! Catalan numbers.

use tertius_integers::Integer;

use crate::binomial::binomial;

/// Computes the nth Catalan number C_n.
///
/// C_n = C(2n, n) / (n + 1) = (2n)! / ((n+1)! * n!)
///
/// Catalan numbers count many combinatorial objects:
/// - Valid sequences of n pairs of parentheses
/// - Full binary trees with n+1 leaves
/// - Paths in a grid from (0,0) to (n,n) not crossing the diagonal
/// - Triangulations of a convex polygon with n+2 sides
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::catalan;
/// use tertius_integers::Integer;
///
/// assert_eq!(catalan(0), Integer::from(1i64));
/// assert_eq!(catalan(1), Integer::from(1i64));
/// assert_eq!(catalan(2), Integer::from(2i64));
/// assert_eq!(catalan(3), Integer::from(5i64));
/// assert_eq!(catalan(4), Integer::from(14i64));
/// assert_eq!(catalan(5), Integer::from(42i64));
/// ```
pub fn catalan(n: u64) -> Integer {
    // C_n = C(2n, n) / (n + 1)
    binomial(2 * n, n) / Integer::from((n + 1) as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_catalan() {
        let expected = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(catalan(n as u64), Integer::from(e as i64));
        }
    }

    #[test]
    fn test_catalan_recurrence() {
        // C_n = sum_{i=0}^{n-1} C_i * C_{n-1-i}
        for n in 1..10 {
            let mut sum = Integer::from(0i64);
            for i in 0..n {
                sum = sum + catalan(i) * catalan(n - 1 - i);
            }
            assert_eq!(catalan(n), sum);
        }
    }

    #[test]
    fn test_catalan_large() {
        // C_20 = 6564120420
        assert_eq!(catalan(20), Integer::from(6564120420i64));
    }
}
