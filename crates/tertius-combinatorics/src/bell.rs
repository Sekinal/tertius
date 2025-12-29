//! Bell numbers.

use num_traits::One;
use tertius_integers::Integer;

/// Computes the nth Bell number B_n.
///
/// B_n is the number of ways to partition a set of n elements.
///
/// B_n = sum_{k=0}^{n} S(n, k) where S(n, k) is the Stirling number of the second kind.
///
/// Alternatively: B_{n+1} = sum_{k=0}^{n} C(n, k) * B_k
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::bell;
/// use tertius_integers::Integer;
///
/// assert_eq!(bell(0), Integer::from(1i64));
/// assert_eq!(bell(1), Integer::from(1i64));
/// assert_eq!(bell(2), Integer::from(2i64));
/// assert_eq!(bell(3), Integer::from(5i64));
/// assert_eq!(bell(4), Integer::from(15i64));
/// assert_eq!(bell(5), Integer::from(52i64));
/// ```
pub fn bell(n: u64) -> Integer {
    if n == 0 {
        return Integer::one();
    }

    // Use the Bell triangle (Aitken's array)
    // Row 0: 1
    // Row 1: 1, 2
    // Row 2: 2, 3, 5
    // Row 3: 5, 7, 10, 15
    // ...
    // B_n is the first element of row n (or last of row n-1)

    let n = n as usize;
    let mut row = vec![Integer::one()]; // Row 0

    for _ in 0..n {
        let mut next_row = Vec::with_capacity(row.len() + 1);
        // First element of next row is last element of current row
        next_row.push(row.last().unwrap().clone());
        // Each subsequent element is sum of previous element and corresponding element in current row
        for i in 0..row.len() {
            let new_val = next_row[i].clone() + row[i].clone();
            next_row.push(new_val);
        }
        row = next_row;
    }

    // B_n is the last element of row n-1, which after n iterations is row[0]
    row[0].clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bell() {
        let expected = [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(bell(n as u64), Integer::from(e as i64));
        }
    }

    #[test]
    fn test_bell_sum_of_stirling() {
        use crate::stirling::stirling2;

        for n in 0..10 {
            let mut sum = Integer::from(0i64);
            for k in 0..=n {
                sum = sum + stirling2(n, k);
            }
            assert_eq!(bell(n), sum);
        }
    }
}
