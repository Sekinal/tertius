//! Stirling numbers of the first and second kind.

use num_traits::Zero;
use tertius_integers::Integer;

/// Computes the (signed) Stirling number of the first kind s(n, k).
///
/// s(n, k) counts the number of permutations of n elements with exactly k cycles,
/// with sign (-1)^(n-k).
///
/// Recurrence: s(n, k) = s(n-1, k-1) - (n-1) * s(n-1, k)
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::stirling1;
/// use tertius_integers::Integer;
///
/// // s(4, 2) = 11 (from OEIS A008275)
/// assert_eq!(stirling1(4, 2), Integer::from(11i64));
///
/// // s(n, n) = 1
/// assert_eq!(stirling1(5, 5), Integer::from(1i64));
/// ```
pub fn stirling1(n: u64, k: u64) -> Integer {
    // Base cases
    if n == 0 && k == 0 {
        return Integer::from(1i64);
    }
    if n == 0 || k == 0 || k > n {
        return Integer::zero();
    }
    if n == k {
        return Integer::from(1i64);
    }

    // Build table using dynamic programming
    let n = n as usize;
    let k = k as usize;

    let mut prev = vec![Integer::zero(); k + 1];
    let mut curr = vec![Integer::zero(); k + 1];

    prev[0] = Integer::from(1i64); // s(0, 0) = 1

    for i in 1..=n {
        curr[0] = Integer::zero();
        for j in 1..=k.min(i) {
            // s(i, j) = s(i-1, j-1) - (i-1) * s(i-1, j)
            let term1 = prev[j - 1].clone();
            let term2 = Integer::from((i - 1) as i64) * prev[j].clone();
            curr[j] = term1 - term2;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[k].clone()
}

/// Computes the Stirling number of the second kind S(n, k).
///
/// S(n, k) counts the number of ways to partition a set of n elements
/// into exactly k non-empty subsets.
///
/// Recurrence: S(n, k) = k * S(n-1, k) + S(n-1, k-1)
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::stirling2;
/// use tertius_integers::Integer;
///
/// // S(4, 2) = 7
/// assert_eq!(stirling2(4, 2), Integer::from(7i64));
///
/// // S(n, 1) = 1 for n >= 1
/// assert_eq!(stirling2(5, 1), Integer::from(1i64));
///
/// // S(n, n) = 1
/// assert_eq!(stirling2(5, 5), Integer::from(1i64));
/// ```
pub fn stirling2(n: u64, k: u64) -> Integer {
    // Base cases
    if n == 0 && k == 0 {
        return Integer::from(1i64);
    }
    if n == 0 || k == 0 || k > n {
        return Integer::zero();
    }
    if k == 1 || n == k {
        return Integer::from(1i64);
    }

    // Build table using dynamic programming
    let n = n as usize;
    let k = k as usize;

    let mut prev = vec![Integer::zero(); k + 1];
    let mut curr = vec![Integer::zero(); k + 1];

    prev[0] = Integer::from(1i64); // S(0, 0) = 1

    for i in 1..=n {
        curr[0] = Integer::zero();
        for j in 1..=k.min(i) {
            // S(i, j) = j * S(i-1, j) + S(i-1, j-1)
            let term1 = Integer::from(j as i64) * prev[j].clone();
            let term2 = prev[j - 1].clone();
            curr[j] = term1 + term2;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[k].clone()
}

/// Computes the unsigned Stirling number of the first kind |s(n, k)|.
///
/// This counts the number of permutations of n elements with exactly k cycles.
///
/// # Examples
///
/// ```ignore
/// // Not currently exported publicly
/// use tertius_combinatorics::stirling::unsigned_stirling1;
/// use tertius_integers::Integer;
///
/// // |s(4, 2)| = 11
/// assert_eq!(unsigned_stirling1(4, 2), Integer::from(11i64));
/// ```
pub fn unsigned_stirling1(n: u64, k: u64) -> Integer {
    stirling1(n, k).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stirling1() {
        // s(n, n) = 1
        for n in 0..10 {
            assert_eq!(stirling1(n, n), Integer::from(1i64));
        }

        // s(n, 0) = 0 for n > 0
        for n in 1..10 {
            assert_eq!(stirling1(n, 0), Integer::zero());
        }

        // s(0, 0) = 1
        assert_eq!(stirling1(0, 0), Integer::from(1i64));

        // Known values from OEIS A008275
        assert_eq!(stirling1(3, 1), Integer::from(2i64));
        assert_eq!(stirling1(3, 2), Integer::from(-3i64));
        assert_eq!(stirling1(4, 2), Integer::from(11i64));
        assert_eq!(stirling1(4, 3), Integer::from(-6i64));
    }

    #[test]
    fn test_stirling2() {
        // S(n, 1) = 1 for n >= 1
        for n in 1..10 {
            assert_eq!(stirling2(n, 1), Integer::from(1i64));
        }

        // S(n, n) = 1
        for n in 0..10 {
            assert_eq!(stirling2(n, n), Integer::from(1i64));
        }

        // S(n, 2) = 2^(n-1) - 1 for n >= 2
        assert_eq!(stirling2(3, 2), Integer::from(3i64));
        assert_eq!(stirling2(4, 2), Integer::from(7i64));
        assert_eq!(stirling2(5, 2), Integer::from(15i64));

        // S(0, 0) = 1
        assert_eq!(stirling2(0, 0), Integer::from(1i64));

        // Known values
        assert_eq!(stirling2(4, 3), Integer::from(6i64));
        assert_eq!(stirling2(5, 3), Integer::from(25i64));
    }

    #[test]
    fn test_stirling_relation() {
        // x^n = sum_k S(n,k) * x(x-1)...(x-k+1)
        // Testing: sum_k S(n,k) * (falling factorial evaluated at some x)
        // This is complex, so we just verify the recurrence
        for n in 2..8 {
            for k in 1..n {
                let lhs = stirling2(n, k);
                let rhs = Integer::from(k as i64) * stirling2(n - 1, k) + stirling2(n - 1, k - 1);
                assert_eq!(lhs, rhs);
            }
        }
    }
}
