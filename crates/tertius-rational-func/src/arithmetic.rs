//! Arithmetic operations for rational functions.
//!
//! Implements field operations: addition, subtraction, multiplication, division.

use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::RationalFunction;
use tertius_rings::traits::Field;

impl<K: Field> Add for RationalFunction<K> {
    type Output = Self;

    /// Adds two rational functions.
    ///
    /// a/b + c/d = (ad + bc) / bd
    fn add(self, other: Self) -> Self::Output {
        self.add_ref(&other)
    }
}

impl<K: Field> Add<&RationalFunction<K>> for RationalFunction<K> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        self.add_ref(other)
    }
}

impl<K: Field> RationalFunction<K> {
    /// Adds two rational functions by reference.
    pub fn add_ref(&self, other: &Self) -> Self {
        // a/b + c/d = (ad + bc) / bd
        let num = self
            .numerator()
            .mul(other.denominator())
            .add(&other.numerator().mul(self.denominator()));
        let den = self.denominator().mul(other.denominator());

        Self::new(num, den)
    }

    /// Subtracts another rational function from this one.
    pub fn sub_ref(&self, other: &Self) -> Self {
        self.add_ref(&other.neg())
    }

    /// Multiplies two rational functions.
    pub fn mul_ref(&self, other: &Self) -> Self {
        let num = self.numerator().mul(other.numerator());
        let den = self.denominator().mul(other.denominator());

        Self::new(num, den)
    }

    /// Divides this rational function by another.
    ///
    /// # Panics
    ///
    /// Panics if `other` is zero.
    pub fn div_ref(&self, other: &Self) -> Self {
        if other.is_zero() {
            panic!("division by zero");
        }

        let num = self.numerator().mul(other.denominator());
        let den = self.denominator().mul(other.numerator());

        Self::new(num, den)
    }

    /// Computes the multiplicative inverse.
    ///
    /// Returns `None` if this is zero.
    pub fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Self::new(
                self.denominator().clone(),
                self.numerator().clone(),
            ))
        }
    }

    /// Raises to a non-negative integer power.
    pub fn pow(&self, n: u32) -> Self {
        if n == 0 {
            return Self::one();
        }

        let mut result = Self::one();
        let mut base = self.clone();
        let mut exp = n;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul_ref(&base);
            }
            base = base.mul_ref(&base);
            exp >>= 1;
        }

        result
    }
}

impl<K: Field> Sub for RationalFunction<K> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        self.sub_ref(&other)
    }
}

impl<K: Field> Sub<&RationalFunction<K>> for RationalFunction<K> {
    type Output = Self;

    fn sub(self, other: &Self) -> Self::Output {
        self.sub_ref(other)
    }
}

impl<K: Field> Mul for RationalFunction<K> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        self.mul_ref(&other)
    }
}

impl<K: Field> Mul<&RationalFunction<K>> for RationalFunction<K> {
    type Output = Self;

    fn mul(self, other: &Self) -> Self::Output {
        self.mul_ref(other)
    }
}

impl<K: Field> Div for RationalFunction<K> {
    type Output = Self;

    fn div(self, other: Self) -> Self::Output {
        self.div_ref(&other)
    }
}

impl<K: Field> Div<&RationalFunction<K>> for RationalFunction<K> {
    type Output = Self;

    fn div(self, other: &Self) -> Self::Output {
        self.div_ref(other)
    }
}

impl<K: Field> Neg for RationalFunction<K> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        RationalFunction::neg(&self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_poly::dense::DensePoly;
    use tertius_rings::rationals::Q;
    use tertius_rings::traits::Ring;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    fn rf(num: &[i64], den: &[i64]) -> RationalFunction<Q> {
        RationalFunction::new(poly(num), poly(den))
    }

    #[test]
    fn test_add() {
        // 1/x + 1/x = 2/x
        let a = rf(&[1], &[0, 1]);
        let b = rf(&[1], &[0, 1]);
        let sum = a + b;

        assert_eq!(sum.numerator().coeff(0), q(2));
        assert_eq!(sum.denominator().degree(), 1);
    }

    #[test]
    fn test_add_different_denominators() {
        // 1/x + 1/(x+1) = (2x + 1) / (x(x+1))
        let a = rf(&[1], &[0, 1]); // 1/x
        let b = rf(&[1], &[1, 1]); // 1/(x+1)
        let sum = a + b;

        // Numerator: (x+1) + x = 2x + 1
        // Denominator: x(x+1) = x^2 + x
        assert_eq!(sum.numerator().degree(), 1);
        assert_eq!(sum.numerator().coeff(0), q(1));
        assert_eq!(sum.numerator().coeff(1), q(2));
    }

    #[test]
    fn test_sub() {
        // 1/x - 1/x = 0
        let a = rf(&[1], &[0, 1]);
        let b = rf(&[1], &[0, 1]);
        let diff = a - b;

        assert!(diff.is_zero());
    }

    #[test]
    fn test_mul() {
        // (1/x) * x = 1
        let a = rf(&[1], &[0, 1]); // 1/x
        let b = RationalFunction::from_poly(poly(&[0, 1])); // x
        let prod = a * b;

        assert!(prod.is_polynomial());
        assert!(prod.numerator().leading_coeff().is_one());
        assert_eq!(prod.numerator().degree(), 0);
    }

    #[test]
    fn test_div() {
        // (1/x) / (1/x) = 1
        let a = rf(&[1], &[0, 1]);
        let b = rf(&[1], &[0, 1]);
        let quot = a / b;

        assert!(quot.is_polynomial());
        assert!(quot.numerator().leading_coeff().is_one());
    }

    #[test]
    fn test_inv() {
        // inv(1/x) = x
        let a = rf(&[1], &[0, 1]);
        let a_inv = a.inv().unwrap();

        assert!(a_inv.is_polynomial());
        assert_eq!(a_inv.numerator().degree(), 1);
    }

    #[test]
    fn test_pow() {
        // (1/x)^2 = 1/x^2
        let a = rf(&[1], &[0, 1]);
        let a_sq = a.pow(2);

        assert_eq!(a_sq.numerator().degree(), 0);
        assert_eq!(a_sq.denominator().degree(), 2);
    }

    #[test]
    fn test_field_identity() {
        // a * inv(a) = 1
        let a = rf(&[1, 2], &[3, 0, 1]); // (1 + 2x) / (3 + x^2)
        let a_inv = a.clone().inv().unwrap();
        let prod = a * a_inv;

        assert!(prod.is_polynomial());
        assert!(prod.numerator().leading_coeff().is_one());
        assert_eq!(prod.numerator().degree(), 0);
    }
}
