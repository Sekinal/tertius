//! Factorial and related functions.

use num_traits::One;
use tertius_integers::Integer;

/// Computes n! (n factorial).
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::factorial;
/// use tertius_integers::Integer;
///
/// assert_eq!(factorial(0), Integer::from(1i64));
/// assert_eq!(factorial(5), Integer::from(120i64));
/// assert_eq!(factorial(10), Integer::from(3628800i64));
/// ```
pub fn factorial(n: u64) -> Integer {
    if n <= 1 {
        return Integer::one();
    }

    let mut result = Integer::one();
    for i in 2..=n {
        result = result * Integer::from(i as i64);
    }
    result
}

/// Computes the falling factorial (n)_k = n * (n-1) * ... * (n-k+1).
///
/// Also known as the Pochhammer symbol (with falling index).
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::falling_factorial;
/// use tertius_integers::Integer;
///
/// // (5)_3 = 5 * 4 * 3 = 60
/// assert_eq!(falling_factorial(5, 3), Integer::from(60i64));
///
/// // (n)_0 = 1
/// assert_eq!(falling_factorial(10, 0), Integer::from(1i64));
/// ```
pub fn falling_factorial(n: i64, k: u64) -> Integer {
    if k == 0 {
        return Integer::one();
    }

    let mut result = Integer::one();
    for i in 0..k {
        result = result * Integer::from(n - i as i64);
    }
    result
}

/// Computes the rising factorial n^(k) = n * (n+1) * ... * (n+k-1).
///
/// Also known as the Pochhammer symbol.
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::rising_factorial;
/// use tertius_integers::Integer;
///
/// // 3^(4) = 3 * 4 * 5 * 6 = 360
/// assert_eq!(rising_factorial(3, 4), Integer::from(360i64));
///
/// // n^(0) = 1
/// assert_eq!(rising_factorial(5, 0), Integer::from(1i64));
/// ```
pub fn rising_factorial(n: i64, k: u64) -> Integer {
    if k == 0 {
        return Integer::one();
    }

    let mut result = Integer::one();
    for i in 0..k {
        result = result * Integer::from(n + i as i64);
    }
    result
}

/// Computes the double factorial n!!.
///
/// n!! = n * (n-2) * (n-4) * ... * (2 or 1)
///
/// # Examples
///
/// ```
/// use tertius_combinatorics::double_factorial;
/// use tertius_integers::Integer;
///
/// // 5!! = 5 * 3 * 1 = 15
/// assert_eq!(double_factorial(5), Integer::from(15i64));
///
/// // 6!! = 6 * 4 * 2 = 48
/// assert_eq!(double_factorial(6), Integer::from(48i64));
/// ```
pub fn double_factorial(n: u64) -> Integer {
    if n <= 1 {
        return Integer::one();
    }

    let mut result = Integer::one();
    let mut i = n;
    while i > 1 {
        result = result * Integer::from(i as i64);
        i -= 2;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), Integer::from(1i64));
        assert_eq!(factorial(1), Integer::from(1i64));
        assert_eq!(factorial(5), Integer::from(120i64));
        assert_eq!(factorial(10), Integer::from(3628800i64));
    }

    #[test]
    fn test_falling_factorial() {
        assert_eq!(falling_factorial(5, 0), Integer::from(1i64));
        assert_eq!(falling_factorial(5, 1), Integer::from(5i64));
        assert_eq!(falling_factorial(5, 3), Integer::from(60i64));
        assert_eq!(falling_factorial(5, 5), Integer::from(120i64));
    }

    #[test]
    fn test_rising_factorial() {
        assert_eq!(rising_factorial(1, 0), Integer::from(1i64));
        assert_eq!(rising_factorial(1, 5), Integer::from(120i64)); // 1*2*3*4*5
        assert_eq!(rising_factorial(3, 4), Integer::from(360i64)); // 3*4*5*6
    }

    #[test]
    fn test_double_factorial() {
        assert_eq!(double_factorial(0), Integer::from(1i64));
        assert_eq!(double_factorial(1), Integer::from(1i64));
        assert_eq!(double_factorial(5), Integer::from(15i64));
        assert_eq!(double_factorial(6), Integer::from(48i64));
        assert_eq!(double_factorial(7), Integer::from(105i64)); // 7*5*3*1
    }
}
