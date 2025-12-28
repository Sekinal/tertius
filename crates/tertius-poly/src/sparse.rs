//! Sparse multivariate polynomials.
//!
//! This module provides sparse polynomial representation for
//! efficient handling of polynomials with few non-zero terms.

use tertius_rings::traits::Ring;

use crate::monomial::PackedMonomial;
use crate::ordering::MonomialOrder;

/// A sparse multivariate polynomial.
///
/// Terms are stored as (monomial, coefficient) pairs, sorted by
/// the monomial ordering.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SparsePoly<R: Ring> {
    /// Terms in sorted order (by monomial ordering).
    terms: Vec<(PackedMonomial, R)>,
    /// Number of variables.
    num_vars: usize,
    /// Monomial ordering used for sorting.
    order: MonomialOrder,
}

impl<R: Ring> SparsePoly<R> {
    /// Creates a new polynomial from terms.
    ///
    /// Terms are automatically sorted and combined.
    #[must_use]
    pub fn new(terms: Vec<(PackedMonomial, R)>, num_vars: usize, order: MonomialOrder) -> Self {
        let mut poly = Self {
            terms,
            num_vars,
            order,
        };
        poly.normalize();
        poly
    }

    /// Creates the zero polynomial.
    #[must_use]
    pub fn zero(num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: Vec::new(),
            num_vars,
            order,
        }
    }

    /// Creates the constant polynomial 1.
    #[must_use]
    pub fn one(num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: vec![(PackedMonomial::one(num_vars), R::one())],
            num_vars,
            order,
        }
    }

    /// Creates a constant polynomial.
    #[must_use]
    pub fn constant(c: R, num_vars: usize, order: MonomialOrder) -> Self {
        if c.is_zero() {
            Self::zero(num_vars, order)
        } else {
            Self {
                terms: vec![(PackedMonomial::one(num_vars), c)],
                num_vars,
                order,
            }
        }
    }

    /// Creates a single variable x_i.
    #[must_use]
    pub fn var(i: usize, num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: vec![(PackedMonomial::var(i, num_vars), R::one())],
            num_vars,
            order,
        }
    }

    /// Returns true if this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of terms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Returns true if there are no terms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of variables.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the monomial ordering.
    #[must_use]
    pub fn order(&self) -> MonomialOrder {
        self.order
    }

    /// Returns the terms.
    #[must_use]
    pub fn terms(&self) -> &[(PackedMonomial, R)] {
        &self.terms
    }

    /// Returns the leading monomial.
    #[must_use]
    pub fn leading_monomial(&self) -> Option<&PackedMonomial> {
        self.terms.first().map(|(m, _)| m)
    }

    /// Returns the leading coefficient.
    #[must_use]
    pub fn leading_coeff(&self) -> Option<&R> {
        self.terms.first().map(|(_, c)| c)
    }

    /// Returns the leading term (monomial, coefficient).
    #[must_use]
    pub fn leading_term(&self) -> Option<&(PackedMonomial, R)> {
        self.terms.first()
    }

    /// Sorts terms and combines like terms.
    fn normalize(&mut self) {
        // Sort by monomial order (descending for leading term first)
        self.terms.sort_by(|a, b| {
            self.order.compare(&b.0, &a.0, self.num_vars)
        });

        // Combine like terms
        let mut i = 0;
        while i < self.terms.len() {
            let j = i + 1;
            while j < self.terms.len() && self.terms[i].0 == self.terms[j].0 {
                let c = self.terms.remove(j).1;
                self.terms[i].1 = self.terms[i].1.clone() + c;
            }
            if self.terms[i].1.is_zero() {
                self.terms.remove(i);
            } else {
                i += 1;
            }
        }
    }

    /// Adds two polynomials.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.num_vars, other.num_vars);
        assert!(self.order == other.order);

        let mut terms = self.terms.clone();
        terms.extend(other.terms.clone());

        Self::new(terms, self.num_vars, self.order)
    }

    /// Negates a polynomial.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            terms: self.terms.iter().map(|(m, c)| (*m, -c.clone())).collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Subtracts two polynomials.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiplies two polynomials (schoolbook algorithm).
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.num_vars, other.num_vars);
        assert!(self.order == other.order);

        if self.is_zero() || other.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        let mut terms = Vec::with_capacity(self.len() * other.len());

        for (m1, c1) in &self.terms {
            for (m2, c2) in &other.terms {
                let m = m1.mul(m2);
                let c = c1.clone() * c2.clone();
                terms.push((m, c));
            }
        }

        Self::new(terms, self.num_vars, self.order)
    }

    /// Multiplies by a scalar.
    #[must_use]
    pub fn scale(&self, c: &R) -> Self {
        if c.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        Self {
            terms: self.terms.iter().map(|(m, x)| (*m, x.clone() * c.clone())).collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Multiplies by a monomial.
    #[must_use]
    pub fn mul_monomial(&self, m: &PackedMonomial, c: &R) -> Self {
        if c.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        Self {
            terms: self.terms.iter().map(|(m2, c2)| (m.mul(m2), c2.clone() * c.clone())).collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Computes the total degree.
    #[must_use]
    pub fn total_degree(&self) -> u32 {
        self.terms
            .iter()
            .map(|(m, _)| m.total_degree(self.num_vars))
            .max()
            .unwrap_or(0)
    }

    /// Returns the degree in a specific variable.
    ///
    /// This is the maximum exponent of `var` across all terms.
    #[must_use]
    pub fn degree_in(&self, var: usize) -> u32 {
        self.terms
            .iter()
            .map(|(m, _)| m.exponent(var))
            .max()
            .unwrap_or(0)
    }

    /// Evaluates the polynomial at a specific value for one variable.
    ///
    /// Returns a polynomial in `num_vars - 1` variables where variable `var`
    /// has been substituted with `value`.
    ///
    /// # Arguments
    /// * `var` - The variable index to substitute (0-indexed)
    /// * `value` - The value to substitute
    ///
    /// # Example
    /// For f(x,y,z) = x²y + 3xyz + z, evaluating at y=2 gives:
    /// f(x,2,z) = 2x² + 6xz + z (a polynomial in x and z only)
    #[must_use]
    pub fn partial_eval(&self, var: usize, value: &R) -> Self {
        if self.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        // Group terms by their exponent in the evaluated variable
        // Then compute value^exp * (terms with that exponent, variable removed)
        let mut result_terms: Vec<(PackedMonomial, R)> = Vec::new();

        for (mono, coeff) in &self.terms {
            let exp = mono.exponent(var);

            // Compute value^exp
            let mut power = R::one();
            for _ in 0..exp {
                power = power * value.clone();
            }

            // Create new monomial without the evaluated variable
            // (keep same num_vars for simplicity - the variable just has exponent 0)
            let mut new_exps = mono.exponents(self.num_vars);
            new_exps[var] = 0;
            let new_mono = PackedMonomial::from_exponents(&new_exps);

            // Add coefficient * value^exp
            let new_coeff = coeff.clone() * power;
            result_terms.push((new_mono, new_coeff));
        }

        Self::new(result_terms, self.num_vars, self.order)
    }

    /// Returns the leading coefficient as a polynomial in the remaining variables.
    ///
    /// Given a polynomial f viewed as an element of R[x_1,...,x_{n-1}][x_main],
    /// this returns the leading coefficient (highest degree coefficient in x_main)
    /// as a polynomial in the other variables.
    ///
    /// # Arguments
    /// * `main_var` - The "main" variable index
    ///
    /// # Example
    /// For f(x,y,z) = (y+1)x² + 3yz·x + z where main_var=0 (x is main),
    /// the leading coefficient is (y+1), a polynomial in y and z.
    #[must_use]
    pub fn leading_coeff_poly(&self, main_var: usize) -> Self {
        if self.is_zero() {
            return Self::zero(self.num_vars, self.order);
        }

        let main_degree = self.degree_in(main_var);
        self.coeff_poly(main_var, main_degree)
    }

    /// Extracts the coefficient of x_var^power as a polynomial in the other variables.
    ///
    /// # Arguments
    /// * `var` - The variable index
    /// * `power` - The power to extract the coefficient for
    ///
    /// # Returns
    /// A polynomial in the same number of variables where `var` always has exponent 0.
    #[must_use]
    pub fn coeff_poly(&self, var: usize, power: u32) -> Self {
        let mut result_terms: Vec<(PackedMonomial, R)> = Vec::new();

        for (mono, coeff) in &self.terms {
            if mono.exponent(var) == power {
                // Create monomial with var exponent set to 0
                let mut new_exps = mono.exponents(self.num_vars);
                new_exps[var] = 0;
                let new_mono = PackedMonomial::from_exponents(&new_exps);
                result_terms.push((new_mono, coeff.clone()));
            }
        }

        Self::new(result_terms, self.num_vars, self.order)
    }

    /// Converts a univariate-like sparse polynomial to dense form.
    ///
    /// This is useful when the polynomial is effectively univariate in one variable
    /// (all other variables have exponent 0).
    ///
    /// # Arguments
    /// * `var` - The main variable
    ///
    /// # Returns
    /// A vector of coefficients where index i contains the coefficient of x_var^i.
    /// Returns None if any term has non-zero exponents in other variables.
    #[must_use]
    pub fn to_univariate_coeffs(&self, var: usize) -> Option<Vec<R>> {
        if self.is_zero() {
            return Some(vec![]);
        }

        let degree = self.degree_in(var);
        let mut coeffs = vec![R::zero(); degree as usize + 1];

        for (mono, coeff) in &self.terms {
            // Check that all other variables have exponent 0
            for i in 0..self.num_vars {
                if i != var && mono.exponent(i) != 0 {
                    return None;
                }
            }
            let exp = mono.exponent(var) as usize;
            coeffs[exp] = coeffs[exp].clone() + coeff.clone();
        }

        Some(coeffs)
    }

    /// Creates a sparse polynomial from univariate coefficients.
    ///
    /// # Arguments
    /// * `coeffs` - Coefficients where index i is the coefficient of x_var^i
    /// * `var` - The variable index
    /// * `num_vars` - Total number of variables
    /// * `order` - Monomial ordering
    #[must_use]
    pub fn from_univariate_coeffs(
        coeffs: &[R],
        var: usize,
        num_vars: usize,
        order: MonomialOrder,
    ) -> Self {
        let mut terms = Vec::new();
        for (i, c) in coeffs.iter().enumerate() {
            if !c.is_zero() {
                let mut exps = vec![0u32; num_vars];
                exps[var] = i as u32;
                terms.push((PackedMonomial::from_exponents(&exps), c.clone()));
            }
        }
        Self::new(terms, num_vars, order)
    }
}

impl<R: Ring> std::fmt::Display for SparsePoly<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let terms: Vec<_> = self
            .terms
            .iter()
            .map(|(m, c)| {
                let mon = m.to_string(self.num_vars);
                if mon == "1" {
                    format!("{c:?}")
                } else {
                    format!("{c:?}*{mon}")
                }
            })
            .collect();

        write!(f, "{}", terms.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_sparse_basic() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let y = SparsePoly::<Q>::var(1, 2, order);

        // x + y
        let sum = x.add(&y);
        assert_eq!(sum.len(), 2);
    }

    #[test]
    fn test_sparse_mul() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let one = SparsePoly::constant(Q::from_integer(1), 2, order);

        // (x + 1)^2 = x^2 + 2x + 1
        let xp1 = x.add(&one);
        let sq = xp1.mul(&xp1);
        assert_eq!(sq.len(), 3);
    }

    #[test]
    fn test_degree_in() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 3, order);
        let y = SparsePoly::<Q>::var(1, 3, order);
        let z = SparsePoly::<Q>::var(2, 3, order);

        // f = x²y + xyz + z² (degree: x->2, y->1, z->2)
        let x2 = x.mul(&x);
        let x2y = x2.mul(&y);
        let xyz = x.mul(&y).mul(&z);
        let z2 = z.mul(&z);
        let f = x2y.add(&xyz).add(&z2);

        assert_eq!(f.degree_in(0), 2); // x
        assert_eq!(f.degree_in(1), 1); // y
        assert_eq!(f.degree_in(2), 2); // z
    }

    #[test]
    fn test_partial_eval() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let y = SparsePoly::<Q>::var(1, 2, order);
        let one = SparsePoly::constant(Q::from_integer(1), 2, order);

        // f = xy + y + 1
        let xy = x.mul(&y);
        let f = xy.add(&y).add(&one);

        // Evaluate at y = 2: f(x, 2) = 2x + 2 + 1 = 2x + 3
        let two = Q::from_integer(2);
        let result = f.partial_eval(1, &two);

        // Should have 2 terms: 2x and 3
        assert_eq!(result.len(), 2);
        assert_eq!(result.degree_in(0), 1); // Still degree 1 in x
        assert_eq!(result.degree_in(1), 0); // No y terms
    }

    #[test]
    fn test_leading_coeff_poly() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let y = SparsePoly::<Q>::var(1, 2, order);
        let one = SparsePoly::constant(Q::from_integer(1), 2, order);

        // f = (y + 1)x² + 3x + y
        // Leading coefficient w.r.t. x is (y + 1)
        let yp1 = y.add(&one);
        let x2 = x.mul(&x);
        let yp1_x2 = yp1.mul(&x2);
        let three = SparsePoly::constant(Q::from_integer(3), 2, order);
        let three_x = three.mul(&x);
        let f = yp1_x2.add(&three_x).add(&y);

        let lc = f.leading_coeff_poly(0); // x is main variable

        // lc should be y + 1
        assert_eq!(lc.len(), 2);
        assert_eq!(lc.degree_in(0), 0); // No x
        assert_eq!(lc.degree_in(1), 1); // Degree 1 in y
    }

    #[test]
    fn test_coeff_poly() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let y = SparsePoly::<Q>::var(1, 2, order);

        // f = x²y + 2xy + 3y + 4
        // coeff of x^1 is 2y
        let x2y = x.mul(&x).mul(&y);
        let two_xy = SparsePoly::constant(Q::from_integer(2), 2, order).mul(&x).mul(&y);
        let three_y = SparsePoly::constant(Q::from_integer(3), 2, order).mul(&y);
        let four = SparsePoly::constant(Q::from_integer(4), 2, order);
        let f = x2y.add(&two_xy).add(&three_y).add(&four);

        let coeff_x1 = f.coeff_poly(0, 1);
        assert_eq!(coeff_x1.len(), 1); // Just 2y
        assert_eq!(coeff_x1.degree_in(1), 1);
    }

    #[test]
    fn test_univariate_conversion() {
        let order = MonomialOrder::Grevlex;

        // Create x² + 2x + 3 as a sparse polynomial in 2 variables (but univariate in x)
        let x = SparsePoly::<Q>::var(0, 2, order);
        let x2 = x.mul(&x);
        let two_x = SparsePoly::constant(Q::from_integer(2), 2, order).mul(&x);
        let three = SparsePoly::constant(Q::from_integer(3), 2, order);
        let f = x2.add(&two_x).add(&three);

        // Convert to univariate coefficients
        let coeffs = f.to_univariate_coeffs(0).unwrap();
        assert_eq!(coeffs.len(), 3);
        assert_eq!(coeffs[0], Q::from_integer(3));
        assert_eq!(coeffs[1], Q::from_integer(2));
        assert_eq!(coeffs[2], Q::from_integer(1));

        // Convert back
        let f2 = SparsePoly::from_univariate_coeffs(&coeffs, 0, 2, order);
        assert_eq!(f, f2);
    }

    #[test]
    fn test_univariate_conversion_fails_for_multivariate() {
        let order = MonomialOrder::Grevlex;
        let x = SparsePoly::<Q>::var(0, 2, order);
        let y = SparsePoly::<Q>::var(1, 2, order);

        // f = x + y is not univariate in x
        let f = x.add(&y);
        assert!(f.to_univariate_coeffs(0).is_none());
    }
}
