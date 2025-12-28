//! Dense univariate polynomials.
//!
//! This module provides high-performance polynomial arithmetic
//! using automatic algorithm selection based on degree.

use tertius_rings::traits::Ring;

/// A dense univariate polynomial.
///
/// Coefficients are stored in ascending degree order.
/// Uses automatic algorithm selection for multiplication:
/// - Degree < 32: Schoolbook
/// - Degree 32-1023: Karatsuba
/// - Degree >= 1024: FFT/NTT
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct DensePoly<R: Ring> {
    /// Coefficients in ascending degree order.
    coeffs: Vec<R>,
}

impl<R: Ring> DensePoly<R> {
    /// Creates a new polynomial from coefficients.
    #[must_use]
    pub fn new(mut coeffs: Vec<R>) -> Self {
        // Normalize: remove trailing zeros
        while coeffs.len() > 1 && coeffs.last().map_or(false, |c| c.is_zero()) {
            coeffs.pop();
        }

        if coeffs.is_empty() {
            coeffs.push(R::zero());
        }

        Self { coeffs }
    }

    /// Creates the zero polynomial.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            coeffs: vec![R::zero()],
        }
    }

    /// Creates the constant polynomial 1.
    #[must_use]
    pub fn one() -> Self {
        Self {
            coeffs: vec![R::one()],
        }
    }

    /// Creates a constant polynomial.
    #[must_use]
    pub fn constant(c: R) -> Self {
        Self::new(vec![c])
    }

    /// Creates the polynomial x.
    #[must_use]
    pub fn x() -> Self {
        Self::new(vec![R::zero(), R::one()])
    }

    /// Creates the monomial c * x^n.
    #[must_use]
    pub fn monomial(c: R, n: usize) -> Self {
        let mut coeffs = vec![R::zero(); n + 1];
        coeffs[n] = c;
        Self::new(coeffs)
    }

    /// Returns the degree of the polynomial.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// Returns true if this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0].is_zero()
    }

    /// Returns the leading coefficient.
    #[must_use]
    pub fn leading_coeff(&self) -> &R {
        self.coeffs.last().unwrap()
    }

    /// Returns the coefficient of x^i.
    #[must_use]
    pub fn coeff(&self, i: usize) -> R {
        self.coeffs.get(i).cloned().unwrap_or_else(R::zero)
    }

    /// Returns all coefficients.
    #[must_use]
    pub fn coeffs(&self) -> &[R] {
        &self.coeffs
    }

    /// Returns mutable access to coefficients.
    pub fn coeffs_mut(&mut self) -> &mut Vec<R> {
        &mut self.coeffs
    }

    /// Evaluates the polynomial at a point using Horner's method.
    #[must_use]
    pub fn eval(&self, x: &R) -> R {
        let mut result = R::zero();
        for c in self.coeffs.iter().rev() {
            result = result * x.clone() + c.clone();
        }
        result
    }

    /// Adds two polynomials.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut result = Vec::with_capacity(len);

        for i in 0..len {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            let b = other.coeffs.get(i).cloned().unwrap_or_else(R::zero);
            result.push(a + b);
        }

        Self::new(result)
    }

    /// Negates a polynomial.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self::new(self.coeffs.iter().map(|c| -c.clone()).collect())
    }

    /// Subtracts two polynomials.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiplies two polynomials.
    ///
    /// Automatically selects the best algorithm based on degree.
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }

        let n = self.degree();
        let m = other.degree();
        let max_deg = n.max(m);

        if max_deg < 32 {
            self.mul_schoolbook(other)
        } else if max_deg < 1024 {
            self.mul_karatsuba(other)
        } else {
            // TODO: Implement FFT multiplication
            self.mul_karatsuba(other)
        }
    }

    /// Schoolbook multiplication: O(n²).
    fn mul_schoolbook(&self, other: &Self) -> Self {
        let n = self.coeffs.len();
        let m = other.coeffs.len();
        let mut result = vec![R::zero(); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                result[i + j] = result[i + j].clone() + self.coeffs[i].clone() * other.coeffs[j].clone();
            }
        }

        Self::new(result)
    }

    /// Karatsuba multiplication: O(n^1.58).
    fn mul_karatsuba(&self, other: &Self) -> Self {
        let n = self.coeffs.len();
        let m = other.coeffs.len();

        // Base case: use schoolbook for small polynomials
        if n < 32 || m < 32 {
            return self.mul_schoolbook(other);
        }

        // Make both polynomials the same size (power of 2)
        let size = (n.max(m)).next_power_of_two();
        let half = size / 2;

        let mut a_coeffs = self.coeffs.clone();
        let mut b_coeffs = other.coeffs.clone();
        a_coeffs.resize(size, R::zero());
        b_coeffs.resize(size, R::zero());

        // Split: a = a0 + a1*x^half, b = b0 + b1*x^half
        let a0 = Self::new(a_coeffs[..half].to_vec());
        let a1 = Self::new(a_coeffs[half..].to_vec());
        let b0 = Self::new(b_coeffs[..half].to_vec());
        let b1 = Self::new(b_coeffs[half..].to_vec());

        // Karatsuba: a*b = z2*x^(2*half) + z1*x^half + z0
        // where z0 = a0*b0, z2 = a1*b1, z1 = (a0+a1)*(b0+b1) - z0 - z2
        let z0 = a0.mul_karatsuba(&b0);
        let z2 = a1.mul_karatsuba(&b1);
        let z1_sum = a0.add(&a1).mul_karatsuba(&b0.add(&b1));
        let z1 = z1_sum.sub(&z0).sub(&z2);

        // Combine
        let mut result = vec![R::zero(); 2 * size - 1];

        for (i, c) in z0.coeffs.iter().enumerate() {
            result[i] = c.clone();
        }

        for (i, c) in z1.coeffs.iter().enumerate() {
            result[i + half] = result[i + half].clone() + c.clone();
        }

        for (i, c) in z2.coeffs.iter().enumerate() {
            result[i + 2 * half] = result[i + 2 * half].clone() + c.clone();
        }

        Self::new(result)
    }

    /// Multiplies by a scalar.
    #[must_use]
    pub fn scale(&self, c: &R) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        Self::new(self.coeffs.iter().map(|x| x.clone() * c.clone()).collect())
    }

    /// Computes the formal derivative.
    #[must_use]
    pub fn derivative(&self) -> Self {
        if self.degree() == 0 {
            return Self::zero();
        }

        let mut result = Vec::with_capacity(self.coeffs.len() - 1);
        for (i, c) in self.coeffs.iter().skip(1).enumerate() {
            result.push(c.mul_by_scalar((i + 1) as i64));
        }

        Self::new(result)
    }

    /// Shifts the polynomial by multiplying by x^n.
    #[must_use]
    pub fn shift(&self, n: usize) -> Self {
        if self.is_zero() || n == 0 {
            return self.clone();
        }

        let mut coeffs = vec![R::zero(); n];
        coeffs.extend(self.coeffs.clone());
        Self::new(coeffs)
    }

    /// Raises the polynomial to a non-negative integer power.
    #[must_use]
    pub fn pow(&self, n: u32) -> Self {
        if n == 0 {
            return Self::one();
        }
        if n == 1 {
            return self.clone();
        }

        let mut result = Self::one();
        let mut base = self.clone();
        let mut exp = n;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exp >>= 1;
        }

        result
    }
}

impl<R: Ring> std::fmt::Display for DensePoly<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut terms = Vec::new();
        for (i, c) in self.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }

            let term = match i {
                0 => format!("{c:?}"),
                1 => format!("{c:?}*x"),
                _ => format!("{c:?}*x^{i}"),
            };
            terms.push(term);
        }

        write!(f, "{}", terms.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_basic_ops() {
        let p = DensePoly::new(vec![Q::from_integer(1), Q::from_integer(2)]); // 1 + 2x
        let q = DensePoly::new(vec![Q::from_integer(3), Q::from_integer(4)]); // 3 + 4x

        let sum = p.add(&q);
        assert_eq!(sum.coeff(0), Q::from_integer(4));
        assert_eq!(sum.coeff(1), Q::from_integer(6));
    }

    #[test]
    fn test_mul_schoolbook() {
        let p = DensePoly::new(vec![Q::from_integer(1), Q::from_integer(2)]); // 1 + 2x
        let q = DensePoly::new(vec![Q::from_integer(3), Q::from_integer(4)]); // 3 + 4x

        // (1 + 2x)(3 + 4x) = 3 + 10x + 8x^2
        let prod = p.mul(&q);
        assert_eq!(prod.coeff(0), Q::from_integer(3));
        assert_eq!(prod.coeff(1), Q::from_integer(10));
        assert_eq!(prod.coeff(2), Q::from_integer(8));
    }

    #[test]
    fn test_eval() {
        // p(x) = 1 + 2x + 3x^2
        let p = DensePoly::new(vec![
            Q::from_integer(1),
            Q::from_integer(2),
            Q::from_integer(3),
        ]);

        // p(2) = 1 + 4 + 12 = 17
        assert_eq!(p.eval(&Q::from_integer(2)), Q::from_integer(17));
    }
}
