//! Integer partitions.

use num_traits::One;
use tertius_integers::Integer;

/// Computes the partition function p(n).
///
/// p(n) is the number of ways to write n as a sum of positive integers,
/// where order doesn't matter.
///
/// For example, p(4) = 5 because 4 can be written as:
/// - 4
/// - 3 + 1
/// - 2 + 2
/// - 2 + 1 + 1
/// - 1 + 1 + 1 + 1
///
/// Uses the recurrence based on Euler's pentagonal number theorem:
/// p(n) = sum_{k≠0} (-1)^(k+1) * p(n - k(3k-1)/2)
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::partition_count;
/// use tertius_integers::Integer;
///
/// assert_eq!(partition_count(0), Integer::from(1i64));
/// assert_eq!(partition_count(1), Integer::from(1i64));
/// assert_eq!(partition_count(4), Integer::from(5i64));
/// assert_eq!(partition_count(10), Integer::from(42i64));
/// ```
pub fn partition_count(n: u64) -> Integer {
    if n == 0 {
        return Integer::one();
    }

    let n = n as usize;

    // Build table using dynamic programming
    let mut p = vec![Integer::from(0i64); n + 1];
    p[0] = Integer::one();

    for i in 1..=n {
        // Use Euler's pentagonal number theorem
        // p(n) = sum_{k≠0} (-1)^(k+1) * p(n - g_k)
        // where g_k = k(3k-1)/2 (generalized pentagonal numbers)

        let mut k = 1i64;
        loop {
            // g_k for positive k
            let g1 = (k * (3 * k - 1) / 2) as usize;
            if g1 > i {
                break;
            }
            let sign = if k % 2 == 1 {
                Integer::one()
            } else {
                Integer::from(-1i64)
            };
            p[i] = p[i].clone() + sign * p[i - g1].clone();

            // g_k for negative k (equivalent to using -k)
            let g2 = (k * (3 * k + 1) / 2) as usize;
            if g2 > i {
                k += 1;
                continue;
            }
            let sign = if k % 2 == 1 {
                Integer::one()
            } else {
                Integer::from(-1i64)
            };
            p[i] = p[i].clone() + sign * p[i - g2].clone();

            k += 1;
        }
    }

    p[n].clone()
}

/// Generates all partitions of n.
///
/// Returns a vector of partitions, where each partition is a vector of positive integers
/// in non-increasing order.
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::partitions;
///
/// let parts = partitions(4);
/// assert_eq!(parts.len(), 5);
/// assert!(parts.contains(&vec![4]));
/// assert!(parts.contains(&vec![3, 1]));
/// assert!(parts.contains(&vec![2, 2]));
/// assert!(parts.contains(&vec![2, 1, 1]));
/// assert!(parts.contains(&vec![1, 1, 1, 1]));
/// ```
pub fn partitions(n: u64) -> Vec<Vec<u64>> {
    if n == 0 {
        return vec![vec![]];
    }

    let mut result = Vec::new();
    let mut current = Vec::new();
    generate_partitions(n, n, &mut current, &mut result);
    result
}

/// Helper function to generate partitions recursively.
fn generate_partitions(n: u64, max_part: u64, current: &mut Vec<u64>, result: &mut Vec<Vec<u64>>) {
    if n == 0 {
        result.push(current.clone());
        return;
    }

    // Add parts from min(n, max_part) down to 1
    for part in (1..=n.min(max_part)).rev() {
        current.push(part);
        generate_partitions(n - part, part, current, result);
        current.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_count() {
        let expected = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(partition_count(n as u64), Integer::from(e as i64));
        }
    }

    #[test]
    fn test_partitions() {
        let parts = partitions(4);
        assert_eq!(parts.len(), 5);
        assert!(parts.contains(&vec![4]));
        assert!(parts.contains(&vec![3, 1]));
        assert!(parts.contains(&vec![2, 2]));
        assert!(parts.contains(&vec![2, 1, 1]));
        assert!(parts.contains(&vec![1, 1, 1, 1]));
    }

    #[test]
    fn test_partitions_count_matches() {
        for n in 0..12 {
            assert_eq!(
                partitions(n as u64).len() as u64,
                partition_count(n as u64).to_i64().unwrap_or(i64::MAX) as u64
            );
        }
    }

    #[test]
    fn test_partitions_sum() {
        // Each partition should sum to n
        for n in 0..10 {
            for part in partitions(n) {
                assert_eq!(part.iter().sum::<u64>(), n);
            }
        }
    }

    #[test]
    fn test_partition_count_large() {
        // p(50) = 204226
        assert_eq!(partition_count(50), Integer::from(204226i64));
    }
}
