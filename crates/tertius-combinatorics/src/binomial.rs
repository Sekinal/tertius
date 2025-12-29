//! Binomial coefficients.

use num_traits::One;
use tertius_integers::Integer;

/// Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
///
/// Uses the efficient formula that avoids computing full factorials:
/// C(n, k) = (n * (n-1) * ... * (n-k+1)) / k!
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::binomial;
/// use tertius_integers::Integer;
///
/// assert_eq!(binomial(5, 2), Integer::from(10i64));
/// assert_eq!(binomial(10, 3), Integer::from(120i64));
/// assert_eq!(binomial(20, 10), Integer::from(184756i64));
/// ```
pub fn binomial(n: u64, k: u64) -> Integer {
    // Handle edge cases
    if k > n {
        return Integer::from(0i64);
    }
    if k == 0 || k == n {
        return Integer::one();
    }

    // Use symmetry: C(n, k) = C(n, n-k)
    let k = k.min(n - k);

    // Compute C(n, k) = n * (n-1) * ... * (n-k+1) / (k * (k-1) * ... * 1)
    // We interleave multiplications and divisions to avoid overflow
    let mut result = Integer::one();

    for i in 0..k {
        result = result * Integer::from((n - i) as i64);
        result = result / Integer::from((i + 1) as i64);
    }

    result
}

/// Computes the generalized binomial coefficient C(n, k) for any real n.
///
/// For non-negative integer k:
/// C(n, k) = n * (n-1) * ... * (n-k+1) / k!
///
/// This allows fractional and negative n.
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::binomial;
/// use tertius_integers::Integer;
///
/// // Standard case
/// assert_eq!(binomial(5, 2), Integer::from(10i64));
/// ```
pub fn binomial_signed(n: i64, k: u64) -> Integer {
    if k == 0 {
        return Integer::one();
    }

    // Use falling factorial: C(n, k) = (n)_k / k!
    let mut result = Integer::one();

    for i in 0..k {
        result = result * Integer::from(n - i as i64);
        result = result / Integer::from((i + 1) as i64);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binomial_basic() {
        assert_eq!(binomial(0, 0), Integer::from(1i64));
        assert_eq!(binomial(5, 0), Integer::from(1i64));
        assert_eq!(binomial(5, 5), Integer::from(1i64));
        assert_eq!(binomial(5, 2), Integer::from(10i64));
        assert_eq!(binomial(10, 3), Integer::from(120i64));
    }

    #[test]
    fn test_binomial_symmetry() {
        for n in 0..15 {
            for k in 0..=n {
                assert_eq!(binomial(n, k), binomial(n, n - k));
            }
        }
    }

    #[test]
    fn test_binomial_pascals_rule() {
        // C(n, k) = C(n-1, k-1) + C(n-1, k)
        for n in 1..15 {
            for k in 1..n {
                let left = binomial(n, k);
                let right = binomial(n - 1, k - 1) + binomial(n - 1, k);
                assert_eq!(left, right);
            }
        }
    }

    #[test]
    fn test_binomial_k_greater_than_n() {
        assert_eq!(binomial(5, 6), Integer::from(0i64));
        assert_eq!(binomial(0, 1), Integer::from(0i64));
    }

    #[test]
    fn test_binomial_large() {
        // C(20, 10) = 184756
        assert_eq!(binomial(20, 10), Integer::from(184756i64));
        // C(30, 15) = 155117520
        assert_eq!(binomial(30, 15), Integer::from(155117520i64));
    }

    #[test]
    fn test_binomial_signed_negative() {
        // C(-1, 2) = (-1)(-2) / 2! = 1
        assert_eq!(binomial_signed(-1, 2), Integer::from(1i64));
        // C(-3, 2) = (-3)(-4) / 2! = 6
        assert_eq!(binomial_signed(-3, 2), Integer::from(6i64));
    }
}
