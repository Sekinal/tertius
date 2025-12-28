//! Core rational function type.
//!
//! A rational function is a quotient of two polynomials P(x)/Q(x).
//! We maintain the invariant that the representation is canonical:
//! - The denominator is monic (leading coefficient = 1)
//! - The numerator and denominator are coprime (gcd = 1)
//! - Zero is represented as 0/1

use tertius_poly::algorithms::gcd::{make_monic, poly_div_rem, poly_gcd};
use tertius_poly::dense::DensePoly;
use tertius_rings::traits::Field;

/// A rational function P(x)/Q(x) over a field K.
///
/// # Invariants
///
/// - `denominator` is always monic (leading coefficient = 1)
/// - `numerator` and `denominator` are coprime (gcd = 1)
/// - Zero is represented as `0 / 1`
///
/// # Example
///
/// ```ignore
/// use tertius_rational_func::RationalFunction;
/// use tertius_rings::rationals::Q;
/// use tertius_poly::dense::DensePoly;
///
/// // Create (x + 1) / (x^2 - 1) = 1 / (x - 1) after normalization
/// let num = DensePoly::new(vec![Q::from_integer(1), Q::from_integer(1)]);
/// let den = DensePoly::new(vec![Q::from_integer(-1), Q::from_integer(0), Q::from_integer(1)]);
/// let f = RationalFunction::new(num, den);
/// ```
#[derive(Clone, Debug)]
pub struct RationalFunction<K: Field> {
    numerator: DensePoly<K>,
    denominator: DensePoly<K>,
}

impl<K: Field> RationalFunction<K> {
    /// Creates a new rational function from numerator and denominator.
    ///
    /// The result is automatically normalized to canonical form.
    ///
    /// # Panics
    ///
    /// Panics if the denominator is zero.
    pub fn new(numerator: DensePoly<K>, denominator: DensePoly<K>) -> Self {
        if denominator.is_zero() {
            panic!("denominator cannot be zero");
        }

        let mut rf = Self {
            numerator,
            denominator,
        };
        rf.normalize();
        rf
    }

    /// Creates a rational function from a polynomial (denominator = 1).
    pub fn from_poly(p: DensePoly<K>) -> Self {
        Self {
            numerator: p,
            denominator: DensePoly::one(),
        }
    }

    /// Creates the zero rational function (0/1).
    pub fn zero() -> Self {
        Self {
            numerator: DensePoly::zero(),
            denominator: DensePoly::one(),
        }
    }

    /// Creates the constant rational function 1/1.
    pub fn one() -> Self {
        Self {
            numerator: DensePoly::one(),
            denominator: DensePoly::one(),
        }
    }

    /// Creates a constant rational function c/1.
    pub fn constant(c: K) -> Self {
        Self {
            numerator: DensePoly::constant(c),
            denominator: DensePoly::one(),
        }
    }

    /// Returns the numerator polynomial.
    pub fn numerator(&self) -> &DensePoly<K> {
        &self.numerator
    }

    /// Returns the denominator polynomial.
    pub fn denominator(&self) -> &DensePoly<K> {
        &self.denominator
    }

    /// Returns true if this is the zero rational function.
    pub fn is_zero(&self) -> bool {
        self.numerator.is_zero()
    }

    /// Returns true if this is a polynomial (denominator = 1).
    pub fn is_polynomial(&self) -> bool {
        self.denominator.degree() == 0 && self.denominator.leading_coeff().is_one()
    }

    /// Returns the polynomial if this is a polynomial, None otherwise.
    pub fn as_polynomial(&self) -> Option<&DensePoly<K>> {
        if self.is_polynomial() {
            Some(&self.numerator)
        } else {
            None
        }
    }

    /// Normalizes the rational function to canonical form.
    ///
    /// This ensures:
    /// 1. The GCD of numerator and denominator is 1
    /// 2. The denominator is monic
    fn normalize(&mut self) {
        // Handle zero numerator
        if self.numerator.is_zero() {
            self.denominator = DensePoly::one();
            return;
        }

        // Compute GCD and divide both by it
        let g = poly_gcd(&self.numerator, &self.denominator);
        if !g.is_zero() && g.degree() > 0 {
            let (num, _) = poly_div_rem(&self.numerator, &g);
            let (den, _) = poly_div_rem(&self.denominator, &g);
            self.numerator = num;
            self.denominator = den;
        }

        // Make denominator monic
        let lead = self.denominator.leading_coeff().clone();
        if !lead.is_one() {
            let lead_inv = lead.inv().expect("field element should have inverse");
            self.numerator = self.numerator.scale(&lead_inv);
            self.denominator = make_monic(&self.denominator);
        }
    }

    /// Returns the degree of the numerator minus degree of denominator.
    ///
    /// This is the "degree at infinity" of the rational function.
    pub fn degree_at_infinity(&self) -> i64 {
        self.numerator.degree() as i64 - self.denominator.degree() as i64
    }

    /// Evaluates the rational function at a point.
    ///
    /// Returns None if the denominator evaluates to zero at that point.
    pub fn eval(&self, x: &K) -> Option<K> {
        let num_val = self.numerator.eval(x);
        let den_val = self.denominator.eval(x);

        if den_val.is_zero() {
            None // pole at x
        } else {
            Some(num_val * den_val.inv()?)
        }
    }

    /// Computes the derivative of the rational function.
    ///
    /// Using the quotient rule: (P/Q)' = (P'Q - PQ') / Q²
    pub fn derivative(&self) -> Self {
        let p = &self.numerator;
        let q = &self.denominator;
        let p_prime = p.derivative();
        let q_prime = q.derivative();

        // P'Q - PQ'
        let num = p_prime.mul(q).sub(&p.mul(&q_prime));
        // Q²
        let den = q.mul(q);

        Self::new(num, den)
    }

    /// Decomposes into polynomial part and proper fraction.
    ///
    /// Returns (poly, proper) where:
    /// - `poly` is the polynomial part (possibly zero)
    /// - `proper` has deg(numerator) < deg(denominator)
    pub fn decompose_proper(&self) -> (DensePoly<K>, Self) {
        if self.numerator.degree() < self.denominator.degree() {
            // Already proper
            return (DensePoly::zero(), self.clone());
        }

        let (q, r) = poly_div_rem(&self.numerator, &self.denominator);

        (q, Self::new(r, self.denominator.clone()))
    }

    /// Negates the rational function.
    pub fn neg(&self) -> Self {
        Self {
            numerator: self.numerator.neg(),
            denominator: self.denominator.clone(),
        }
    }
}

impl<K: Field> PartialEq for RationalFunction<K> {
    fn eq(&self, other: &Self) -> bool {
        // Since both are normalized, we can compare directly
        self.numerator == other.numerator && self.denominator == other.denominator
    }
}

impl<K: Field> Eq for RationalFunction<K> {}

impl<K: Field> std::fmt::Display for RationalFunction<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_polynomial() {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "({}) / ({})", self.numerator, self.denominator)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_basic_creation() {
        // Create (x + 1) / (x^2)
        let num = poly(&[1, 1]); // 1 + x
        let den = poly(&[0, 0, 1]); // x^2
        let rf = RationalFunction::new(num, den);

        assert!(!rf.is_zero());
        assert!(!rf.is_polynomial());
    }

    #[test]
    fn test_normalization_common_factor() {
        // Create (x^2 - 1) / (x - 1) = x + 1 after normalization
        let num = poly(&[-1, 0, 1]); // x^2 - 1 = (x-1)(x+1)
        let den = poly(&[-1, 1]); // x - 1
        let rf = RationalFunction::new(num, den);

        // Should simplify to polynomial x + 1
        assert!(rf.is_polynomial());
        assert_eq!(rf.numerator().degree(), 1);
    }

    #[test]
    fn test_decompose_proper() {
        // (x^3 + 1) / (x + 1) = x^2 - x + 1 + 0/(x+1)
        let num = poly(&[1, 0, 0, 1]); // x^3 + 1
        let den = poly(&[1, 1]); // x + 1
        let rf = RationalFunction::new(num, den);

        let (poly_part, proper) = rf.decompose_proper();

        // x^3 + 1 = (x+1)(x^2 - x + 1)
        // So polynomial part should be x^2 - x + 1 and proper part should be 0
        assert_eq!(poly_part.degree(), 2);
        assert!(proper.is_zero());
    }

    #[test]
    fn test_derivative() {
        // d/dx (1/x) = -1/x^2
        let num = poly(&[1]); // 1
        let den = poly(&[0, 1]); // x
        let rf = RationalFunction::new(num, den);

        let rf_prime = rf.derivative();

        // Should be -1 / x^2
        assert_eq!(rf_prime.numerator().degree(), 0);
        assert_eq!(rf_prime.denominator().degree(), 2);
        assert_eq!(rf_prime.numerator().coeff(0), q(-1));
    }

    #[test]
    #[should_panic(expected = "denominator cannot be zero")]
    fn test_zero_denominator_panics() {
        let num = poly(&[1]);
        let den = DensePoly::zero();
        let _ = RationalFunction::new(num, den);
    }
}
