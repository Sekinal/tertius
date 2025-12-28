//! Laurent series with negative powers.
//!
//! A Laurent series is a generalization of a power series that allows
//! negative powers of x:
//!
//! f(x) = Σᵢ₌₋∞^∞ aᵢxⁱ = ... + a₋₂x⁻² + a₋₁x⁻¹ + a₀ + a₁x + a₂x² + ...
//!
//! In practice, we truncate both ends:
//! f(x) = Σᵢ₌₋ₘ^n aᵢxⁱ

use num_traits::{One, Zero};
use std::collections::BTreeMap;
use tertius_rings::traits::{Field, Ring};

/// A Laurent series with finitely many negative powers.
#[derive(Clone, Debug)]
pub struct LaurentSeries<R: Ring + Clone> {
    /// Coefficients indexed by exponent (can be negative).
    coeffs: BTreeMap<i32, R>,
    /// Minimum exponent (most negative power).
    min_exp: i32,
    /// Maximum exponent (highest positive power).
    max_exp: i32,
}

impl<R: Ring + Clone> LaurentSeries<R> {
    /// Creates a new Laurent series from coefficients.
    ///
    /// Coefficients are given as (exponent, coefficient) pairs.
    pub fn new(terms: Vec<(i32, R)>) -> Self {
        let mut coeffs = BTreeMap::new();
        let mut min_exp = i32::MAX;
        let mut max_exp = i32::MIN;

        for (exp, coeff) in terms {
            if !coeff.is_zero() {
                coeffs.insert(exp, coeff);
                min_exp = min_exp.min(exp);
                max_exp = max_exp.max(exp);
            }
        }

        if coeffs.is_empty() {
            min_exp = 0;
            max_exp = 0;
        }

        Self {
            coeffs,
            min_exp,
            max_exp,
        }
    }

    /// Creates the zero series.
    pub fn zero() -> Self {
        Self::new(vec![])
    }

    /// Creates a constant series.
    pub fn constant(c: R) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self::new(vec![(0, c)])
        }
    }

    /// Creates the series x^n for any integer n.
    pub fn monomial(n: i32) -> Self {
        Self::new(vec![(n, R::one())])
    }

    /// Returns the coefficient of x^n.
    pub fn coeff(&self, n: i32) -> R {
        self.coeffs.get(&n).cloned().unwrap_or_else(R::zero)
    }

    /// Returns the minimum exponent (order of the pole).
    pub fn min_exponent(&self) -> i32 {
        self.min_exp
    }

    /// Returns the maximum exponent.
    pub fn max_exponent(&self) -> i32 {
        self.max_exp
    }

    /// Returns the order (index of first non-zero coefficient).
    pub fn order(&self) -> Option<i32> {
        self.coeffs.keys().next().copied()
    }

    /// Returns true if the series is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Returns true if the series has a pole (negative exponents).
    pub fn has_pole(&self) -> bool {
        self.min_exp < 0
    }

    /// Returns the order of the pole (0 if no pole).
    pub fn pole_order(&self) -> u32 {
        if self.min_exp < 0 {
            (-self.min_exp) as u32
        } else {
            0
        }
    }

    /// Returns all non-zero terms as (exponent, coefficient) pairs.
    pub fn terms(&self) -> impl Iterator<Item = (i32, &R)> {
        self.coeffs.iter().map(|(&e, c)| (e, c))
    }

    /// Returns the principal part (negative powers only).
    pub fn principal_part(&self) -> Self {
        let terms: Vec<_> = self
            .coeffs
            .iter()
            .filter(|(&e, _)| e < 0)
            .map(|(&e, c)| (e, c.clone()))
            .collect();
        Self::new(terms)
    }

    /// Returns the analytic part (non-negative powers only).
    pub fn analytic_part(&self) -> Self {
        let terms: Vec<_> = self
            .coeffs
            .iter()
            .filter(|(&e, _)| e >= 0)
            .map(|(&e, c)| (e, c.clone()))
            .collect();
        Self::new(terms)
    }

    /// Adds two Laurent series.
    pub fn add(&self, other: &Self) -> Self {
        let mut terms = Vec::new();

        let min = self.min_exp.min(other.min_exp);
        let max = self.max_exp.max(other.max_exp);

        for i in min..=max {
            let c = self.coeff(i) + other.coeff(i);
            if !c.is_zero() {
                terms.push((i, c));
            }
        }

        Self::new(terms)
    }

    /// Subtracts two Laurent series.
    pub fn sub(&self, other: &Self) -> Self {
        let mut terms = Vec::new();

        let min = self.min_exp.min(other.min_exp);
        let max = self.max_exp.max(other.max_exp);

        for i in min..=max {
            let c = self.coeff(i) - other.coeff(i);
            if !c.is_zero() {
                terms.push((i, c));
            }
        }

        Self::new(terms)
    }

    /// Scales by a constant.
    pub fn scale(&self, c: &R) -> Self {
        if c.is_zero() {
            return Self::zero();
        }

        let terms: Vec<_> = self
            .coeffs
            .iter()
            .map(|(&e, coeff)| (e, coeff.clone() * c.clone()))
            .filter(|(_, c)| !c.is_zero())
            .collect();

        Self::new(terms)
    }

    /// Multiplies by x^n (shift exponents).
    pub fn shift(&self, n: i32) -> Self {
        let terms: Vec<_> = self
            .coeffs
            .iter()
            .map(|(&e, c)| (e + n, c.clone()))
            .collect();
        Self::new(terms)
    }

    /// Multiplies two Laurent series.
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }

        let mut result = BTreeMap::new();

        for (&e1, c1) in &self.coeffs {
            for (&e2, c2) in &other.coeffs {
                let exp = e1 + e2;
                let prod = c1.clone() * c2.clone();
                result
                    .entry(exp)
                    .and_modify(|c: &mut R| *c = c.clone() + prod.clone())
                    .or_insert(prod);
            }
        }

        // Filter zeros and convert
        let terms: Vec<_> = result
            .into_iter()
            .filter(|(_, c)| !c.is_zero())
            .collect();

        Self::new(terms)
    }
}

/// Operations requiring a field.
impl<R: Field + Clone> LaurentSeries<R> {
    /// Computes the residue (coefficient of x^{-1}).
    pub fn residue(&self) -> R {
        self.coeff(-1)
    }

    /// Computes the multiplicative inverse.
    ///
    /// Requires the series to have a non-zero leading coefficient.
    pub fn inverse(&self, max_terms: usize) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Get the leading term
        let ord = self.order()?;
        let leading = self.coeff(ord);
        if leading.is_zero() {
            return None;
        }

        let leading_inv = leading.inv()?;

        // Normalize: f(x) = x^ord * (a_0 + a_1 x + ...)
        // 1/f = x^{-ord} * 1/(a_0 + a_1 x + ...)

        // For the normalized series g = 1 + higher terms, compute 1/g
        // 1/g = 1 - (g-1) + (g-1)^2 - ...

        // More directly, use the recurrence:
        // If f * h = 1, then h_n = -(1/a_0) * Σᵢ₌₁ⁿ a_i * h_{n-i}

        let mut inv_coeffs = Vec::new();

        for n in 0..max_terms {
            if n == 0 {
                inv_coeffs.push(leading_inv.clone());
            } else {
                let mut sum = R::zero();
                for i in 1..=n {
                    let a_i = self.coeff(ord + i as i32);
                    if !a_i.is_zero() && n - i < inv_coeffs.len() {
                        sum = sum + a_i * inv_coeffs[n - i].clone();
                    }
                }
                inv_coeffs.push(R::zero() - leading_inv.clone() * sum);
            }
        }

        // Convert to terms with shifted exponents
        let terms: Vec<_> = inv_coeffs
            .into_iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| (-ord + i as i32, c))
            .collect();

        Some(Self::new(terms))
    }

    /// Divides two Laurent series.
    pub fn div(&self, other: &Self, max_terms: usize) -> Option<Self> {
        let inv = other.inverse(max_terms)?;
        Some(self.mul(&inv))
    }
}

/// Formatting for display.
impl<R: Ring + Clone + std::fmt::Display> std::fmt::Display for LaurentSeries<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut first = true;
        for (&exp, coeff) in &self.coeffs {
            if first {
                first = false;
            } else {
                write!(f, " + ")?;
            }

            if exp == 0 {
                write!(f, "{}", coeff)?;
            } else if exp == 1 {
                write!(f, "{}*x", coeff)?;
            } else if exp == -1 {
                write!(f, "{}/x", coeff)?;
            } else if exp > 0 {
                write!(f, "{}*x^{}", coeff, exp)?;
            } else {
                write!(f, "{}/x^{}", coeff, -exp)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64, d: i64) -> Q {
        Q::new(n, d)
    }

    #[test]
    fn test_new() {
        let s = LaurentSeries::new(vec![(-1, q(1, 1)), (0, q(2, 1)), (1, q(3, 1))]);

        assert_eq!(s.coeff(-1), q(1, 1));
        assert_eq!(s.coeff(0), q(2, 1));
        assert_eq!(s.coeff(1), q(3, 1));
        assert_eq!(s.min_exponent(), -1);
        assert_eq!(s.max_exponent(), 1);
    }

    #[test]
    fn test_monomial() {
        let x_inv = LaurentSeries::<Q>::monomial(-1);
        assert_eq!(x_inv.coeff(-1), q(1, 1));
        assert!(x_inv.has_pole());
        assert_eq!(x_inv.pole_order(), 1);
    }

    #[test]
    fn test_add() {
        let a = LaurentSeries::new(vec![(0, q(1, 1)), (1, q(2, 1))]);
        let b = LaurentSeries::new(vec![(0, q(3, 1)), (1, q(4, 1))]);
        let sum = a.add(&b);

        assert_eq!(sum.coeff(0), q(4, 1));
        assert_eq!(sum.coeff(1), q(6, 1));
    }

    #[test]
    fn test_mul() {
        // (1 + x) * (1 - x) = 1 - x²
        let a = LaurentSeries::new(vec![(0, q(1, 1)), (1, q(1, 1))]);
        let b = LaurentSeries::new(vec![(0, q(1, 1)), (1, q(-1, 1))]);
        let prod = a.mul(&b);

        assert_eq!(prod.coeff(0), q(1, 1));
        assert_eq!(prod.coeff(1), q(0, 1));
        assert_eq!(prod.coeff(2), q(-1, 1));
    }

    #[test]
    fn test_shift() {
        let a = LaurentSeries::new(vec![(0, q(1, 1)), (1, q(2, 1))]);
        let shifted = a.shift(-2);

        assert_eq!(shifted.coeff(-2), q(1, 1));
        assert_eq!(shifted.coeff(-1), q(2, 1));
        assert_eq!(shifted.min_exponent(), -2);
    }

    #[test]
    fn test_residue() {
        let s = LaurentSeries::new(vec![(-2, q(1, 1)), (-1, q(3, 1)), (0, q(5, 1))]);
        assert_eq!(s.residue(), q(3, 1));
    }

    #[test]
    fn test_principal_part() {
        let s = LaurentSeries::new(vec![(-2, q(1, 1)), (-1, q(2, 1)), (0, q(3, 1)), (1, q(4, 1))]);
        let pp = s.principal_part();

        assert_eq!(pp.coeff(-2), q(1, 1));
        assert_eq!(pp.coeff(-1), q(2, 1));
        assert_eq!(pp.coeff(0), q(0, 1));
    }

    #[test]
    fn test_inverse() {
        // 1/(1 + x) = 1 - x + x² - x³ + ...
        let f = LaurentSeries::new(vec![(0, q(1, 1)), (1, q(1, 1))]);
        let inv = f.inverse(5).unwrap();

        assert_eq!(inv.coeff(0), q(1, 1));
        assert_eq!(inv.coeff(1), q(-1, 1));
        assert_eq!(inv.coeff(2), q(1, 1));
        assert_eq!(inv.coeff(3), q(-1, 1));
    }

    #[test]
    fn test_inverse_pole() {
        // 1/x has inverse x
        let x_inv = LaurentSeries::<Q>::monomial(-1);
        let inv = x_inv.inverse(5).unwrap();

        assert_eq!(inv.coeff(1), q(1, 1));
        assert_eq!(inv.min_exponent(), 1);
    }
}
