//! Bernoulli numbers.

use num_traits::{One, Zero};
use tertius_integers::Integer;

use crate::binomial::binomial;

/// Computes the nth Bernoulli number B_n as a rational number (numerator, denominator).
///
/// Bernoulli numbers are defined by the generating function:
/// x / (e^x - 1) = sum_{n=0}^∞ B_n * x^n / n!
///
/// Key values:
/// - B_0 = 1
/// - B_1 = -1/2 (or +1/2 in some conventions)
/// - B_2 = 1/6
/// - B_odd>1 = 0
///
/// Uses the recurrence: B_n = -1/(n+1) * sum_{k=0}^{n-1} C(n+1, k) * B_k
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::bernoulli;
/// use tertius_integers::Integer;
///
/// let (num, den) = bernoulli(0);
/// assert_eq!(num, Integer::from(1i64));
/// assert_eq!(den, Integer::from(1i64));
///
/// let (num, den) = bernoulli(1);
/// assert_eq!(num, Integer::from(-1i64));
/// assert_eq!(den, Integer::from(2i64));
///
/// let (num, den) = bernoulli(2);
/// assert_eq!(num, Integer::from(1i64));
/// assert_eq!(den, Integer::from(6i64));
/// ```
pub fn bernoulli(n: u64) -> (Integer, Integer) {
    // B_0 = 1
    if n == 0 {
        return (Integer::one(), Integer::one());
    }

    // B_1 = -1/2
    if n == 1 {
        return (Integer::from(-1i64), Integer::from(2i64));
    }

    // B_odd = 0 for odd n > 1
    if n % 2 == 1 {
        return (Integer::zero(), Integer::one());
    }

    // Compute using the recurrence
    // Store computed Bernoulli numbers as (numerator, denominator)
    let n = n as usize;
    let mut b_nums: Vec<Integer> = vec![Integer::one()];
    let mut b_dens: Vec<Integer> = vec![Integer::one()];

    // B_1 = -1/2
    b_nums.push(Integer::from(-1i64));
    b_dens.push(Integer::from(2i64));

    for m in 2..=n {
        if m % 2 == 1 {
            // B_m = 0 for odd m > 1
            b_nums.push(Integer::zero());
            b_dens.push(Integer::one());
            continue;
        }

        // B_m = -1/(m+1) * sum_{k=0}^{m-1} C(m+1, k) * B_k
        // First compute the sum with a common denominator

        // Find LCM of all denominators
        let mut lcm = Integer::one();
        for k in 0..m {
            lcm = lcm_integers(&lcm, &b_dens[k]);
        }

        // Compute sum with common denominator lcm
        let mut sum_num = Integer::zero();
        for k in 0..m {
            let coeff = binomial((m + 1) as u64, k as u64);
            let scale = lcm.clone() / b_dens[k].clone();
            sum_num = sum_num + coeff * b_nums[k].clone() * scale;
        }

        // B_m = -sum / ((m+1) * lcm)
        let num = -sum_num;
        let den = Integer::from((m + 1) as i64) * lcm;

        // Reduce fraction
        let (reduced_num, reduced_den) = reduce_fraction(num, den);

        b_nums.push(reduced_num);
        b_dens.push(reduced_den);
    }

    (b_nums[n].clone(), b_dens[n].clone())
}

/// Computes GCD of two integers.
fn gcd_integers(a: &Integer, b: &Integer) -> Integer {
    a.gcd(b)
}

/// Computes LCM of two integers.
fn lcm_integers(a: &Integer, b: &Integer) -> Integer {
    if a.is_zero() || b.is_zero() {
        return Integer::zero();
    }
    (a.clone() * b.clone()).abs() / gcd_integers(a, b)
}

/// Reduces a fraction to lowest terms.
fn reduce_fraction(num: Integer, den: Integer) -> (Integer, Integer) {
    if num.is_zero() {
        return (Integer::zero(), Integer::one());
    }

    let g = gcd_integers(&num, &den);
    let mut n = num / g.clone();
    let mut d = den / g;

    // Ensure denominator is positive
    if d.is_negative() {
        n = -n;
        d = -d;
    }

    (n, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bernoulli_basic() {
        // B_0 = 1
        let (n, d) = bernoulli(0);
        assert_eq!(n, Integer::from(1i64));
        assert_eq!(d, Integer::from(1i64));

        // B_1 = -1/2
        let (n, d) = bernoulli(1);
        assert_eq!(n, Integer::from(-1i64));
        assert_eq!(d, Integer::from(2i64));

        // B_2 = 1/6
        let (n, d) = bernoulli(2);
        assert_eq!(n, Integer::from(1i64));
        assert_eq!(d, Integer::from(6i64));
    }

    #[test]
    fn test_bernoulli_odd() {
        // B_n = 0 for odd n > 1
        for n in [3, 5, 7, 9, 11].iter() {
            let (num, _) = bernoulli(*n);
            assert!(num.is_zero());
        }
    }

    #[test]
    fn test_bernoulli_even() {
        // B_4 = -1/30
        let (n, d) = bernoulli(4);
        assert_eq!(n, Integer::from(-1i64));
        assert_eq!(d, Integer::from(30i64));

        // B_6 = 1/42
        let (n, d) = bernoulli(6);
        assert_eq!(n, Integer::from(1i64));
        assert_eq!(d, Integer::from(42i64));

        // B_8 = -1/30
        let (n, d) = bernoulli(8);
        assert_eq!(n, Integer::from(-1i64));
        assert_eq!(d, Integer::from(30i64));

        // B_10 = 5/66
        let (n, d) = bernoulli(10);
        assert_eq!(n, Integer::from(5i64));
        assert_eq!(d, Integer::from(66i64));
    }
}
