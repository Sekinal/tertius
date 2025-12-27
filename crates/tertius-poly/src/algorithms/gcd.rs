//! Polynomial GCD algorithms.
//!
//! This module provides algorithms for computing the greatest common
//! divisor of polynomials.

use tertius_rings::traits::{EuclideanDomain, Field};

use crate::dense::DensePoly;

/// Computes the GCD of two polynomials over a field using the Euclidean algorithm.
pub fn poly_gcd<F: Field>(a: &DensePoly<F>, b: &DensePoly<F>) -> DensePoly<F> {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }

    let mut p = a.clone();
    let mut q = b.clone();

    while !q.is_zero() {
        let (_, r) = poly_div_rem(&p, &q);
        p = q;
        q = r;
    }

    // Make monic
    make_monic(&p)
}

/// Divides polynomial a by b, returning (quotient, remainder).
///
/// Requires the coefficient ring to be a field.
pub fn poly_div_rem<F: Field>(a: &DensePoly<F>, b: &DensePoly<F>) -> (DensePoly<F>, DensePoly<F>) {
    if b.is_zero() {
        panic!("division by zero polynomial");
    }

    if a.degree() < b.degree() {
        return (DensePoly::zero(), a.clone());
    }

    let b_lead_inv = b.leading_coeff().inv().expect("field element should have inverse");
    let mut quotient = vec![F::zero(); a.degree() - b.degree() + 1];
    let mut remainder = a.coeffs().to_vec();

    while remainder.len() >= b.coeffs().len() {
        let deg_diff = remainder.len() - b.coeffs().len();
        let coeff = remainder.last().unwrap().clone() * b_lead_inv.clone();

        quotient[deg_diff] = coeff.clone();

        for (i, bc) in b.coeffs().iter().enumerate() {
            remainder[deg_diff + i] = remainder[deg_diff + i].clone() - coeff.clone() * bc.clone();
        }

        // Remove trailing zeros
        while remainder.len() > 1 && remainder.last().map_or(false, |c| c.is_zero()) {
            remainder.pop();
        }

        if remainder.len() == 1 && remainder[0].is_zero() {
            break;
        }
    }

    (DensePoly::new(quotient), DensePoly::new(remainder))
}

/// Makes a polynomial monic (leading coefficient = 1).
pub fn make_monic<F: Field>(p: &DensePoly<F>) -> DensePoly<F> {
    if p.is_zero() {
        return p.clone();
    }

    let lead_inv = p.leading_coeff().inv().expect("field element should have inverse");
    p.scale(&lead_inv)
}

/// Computes the content of a polynomial (GCD of all coefficients).
pub fn content<R: EuclideanDomain>(p: &DensePoly<R>) -> R {
    if p.is_zero() {
        return R::zero();
    }

    p.coeffs()
        .iter()
        .cloned()
        .reduce(|a, b| a.gcd(&b))
        .unwrap_or_else(R::zero)
}

/// Computes the primitive part of a polynomial (divided by content).
pub fn primitive_part<R: EuclideanDomain>(p: &DensePoly<R>) -> DensePoly<R> {
    let c = content(p);
    if c.is_zero() || c.is_one() {
        return p.clone();
    }

    let coeffs: Vec<R> = p.coeffs()
        .iter()
        .map(|x| x.div(&c))
        .collect();

    DensePoly::new(coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;
    use tertius_rings::traits::Ring;

    #[test]
    fn test_poly_div_rem() {
        // (x^2 + 2x + 1) / (x + 1) = (x + 1), remainder 0
        let a = DensePoly::new(vec![
            Q::from_integer(1),
            Q::from_integer(2),
            Q::from_integer(1),
        ]);
        let b = DensePoly::new(vec![Q::from_integer(1), Q::from_integer(1)]);

        let (q, r) = poly_div_rem(&a, &b);

        assert_eq!(q.coeff(0), Q::from_integer(1));
        assert_eq!(q.coeff(1), Q::from_integer(1));
        assert!(r.is_zero());
    }

    #[test]
    fn test_poly_gcd() {
        // gcd(x^2 - 1, x^2 - 2x + 1) = x - 1
        // x^2 - 1 = (x-1)(x+1)
        // x^2 - 2x + 1 = (x-1)^2
        let a = DensePoly::new(vec![
            Q::from_integer(-1),
            Q::from_integer(0),
            Q::from_integer(1),
        ]);
        let b = DensePoly::new(vec![
            Q::from_integer(1),
            Q::from_integer(-2),
            Q::from_integer(1),
        ]);

        let g = poly_gcd(&a, &b);

        // Should be monic, degree 1
        assert_eq!(g.degree(), 1);
        assert!(g.leading_coeff().is_one());
    }
}
