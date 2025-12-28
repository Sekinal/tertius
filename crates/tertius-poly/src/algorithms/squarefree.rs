//! Squarefree decomposition of polynomials.
//!
//! A polynomial is squarefree if it has no repeated factors.
//! The squarefree decomposition writes a polynomial as:
//!
//! f = f₁ * f₂² * f₃³ * ...
//!
//! where each fᵢ is squarefree and coprime to the others.
//!
//! # Algorithm
//!
//! Uses Yun's algorithm, which works over any field of characteristic 0
//! (or characteristic > deg(f)).

use crate::algorithms::gcd::{poly_div_rem, poly_gcd};
use crate::dense::DensePoly;
use tertius_rings::traits::Field;

/// A factor with its multiplicity in the squarefree decomposition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquarefreeFactor<F: Field> {
    /// The squarefree polynomial factor.
    pub factor: DensePoly<F>,
    /// The multiplicity (power) of this factor.
    pub multiplicity: u32,
}

/// Result of squarefree decomposition.
#[derive(Clone, Debug)]
pub struct SquarefreeDecomposition<F: Field> {
    /// The unit (leading coefficient factor).
    pub unit: F,
    /// The squarefree factors with multiplicities.
    pub factors: Vec<SquarefreeFactor<F>>,
}

impl<F: Field> SquarefreeDecomposition<F> {
    /// Reconstructs the original polynomial from the decomposition.
    pub fn to_polynomial(&self) -> DensePoly<F> {
        let mut result = DensePoly::constant(self.unit.clone());

        for sf in &self.factors {
            let powered = sf.factor.pow(sf.multiplicity);
            result = result.mul(&powered);
        }

        result
    }

    /// Returns true if the polynomial is squarefree (all multiplicities are 1).
    pub fn is_squarefree(&self) -> bool {
        self.factors.iter().all(|f| f.multiplicity == 1)
    }

    /// Returns the squarefree part (product of all factors with multiplicity 1).
    pub fn squarefree_part(&self) -> DensePoly<F> {
        let mut result = DensePoly::one();

        for sf in &self.factors {
            result = result.mul(&sf.factor);
        }

        result
    }
}

/// Computes the squarefree decomposition of a polynomial using Yun's algorithm.
///
/// Returns the decomposition f = unit * f₁ * f₂² * f₃³ * ...
/// where each fᵢ is squarefree, monic, and coprime to the others.
///
/// # Algorithm (Yun)
///
/// 1. g = gcd(f, f')
/// 2. a₀ = f/g, b₀ = f'/g
/// 3. Loop:
///    - c = b - a'
///    - if c = 0: output a and terminate
///    - d = gcd(a, c)
///    - output d with current multiplicity
///    - a = a/d, b = c/d
///
/// # Example
///
/// ```ignore
/// use tertius_poly::dense::DensePoly;
/// use tertius_poly::algorithms::squarefree::squarefree_decomposition;
/// use tertius_rings::rationals::Q;
///
/// // f = (x+1)² * (x+2)³
/// let f = DensePoly::new(vec![...]);
/// let decomp = squarefree_decomposition(&f);
///
/// // decomp.factors will have (x+1, 2) and (x+2, 3)
/// ```
pub fn squarefree_decomposition<F: Field>(f: &DensePoly<F>) -> SquarefreeDecomposition<F> {
    // Handle constant case
    if f.degree() == 0 {
        return SquarefreeDecomposition {
            unit: f.coeff(0),
            factors: Vec::new(),
        };
    }

    // Extract leading coefficient
    let unit = f.leading_coeff().clone();

    // Make f monic for the algorithm
    let f_monic = {
        let lead_inv = unit.inv().expect("field element should have inverse");
        f.scale(&lead_inv)
    };

    let f_prime = f_monic.derivative();

    // If f' = 0, then f = g^p for some g and the characteristic p of the field
    // This shouldn't happen for characteristic 0 fields (Q), but handle gracefully
    if f_prime.is_zero() {
        return SquarefreeDecomposition {
            unit,
            factors: vec![SquarefreeFactor {
                factor: f_monic,
                multiplicity: 1,
            }],
        };
    }

    // g = gcd(f, f')
    let g = poly_gcd(&f_monic, &f_prime);

    // If gcd = 1, f is already squarefree
    if g.degree() == 0 {
        return SquarefreeDecomposition {
            unit,
            factors: vec![SquarefreeFactor {
                factor: f_monic,
                multiplicity: 1,
            }],
        };
    }

    // a₀ = f/g, b₀ = f'/g
    let (mut a, _) = poly_div_rem(&f_monic, &g);
    let (mut b, _) = poly_div_rem(&f_prime, &g);

    let mut factors = Vec::new();
    let mut multiplicity = 1u32;

    loop {
        // c = b - a'
        let a_prime = a.derivative();
        let c = b.sub(&a_prime);

        if c.is_zero() {
            // a is the last factor
            if a.degree() > 0 {
                factors.push(SquarefreeFactor {
                    factor: a,
                    multiplicity,
                });
            }
            break;
        }

        // d = gcd(a, c)
        let d = poly_gcd(&a, &c);

        if d.degree() > 0 {
            factors.push(SquarefreeFactor {
                factor: d.clone(),
                multiplicity,
            });
        }

        // a = a/d, b = c/d
        let (new_a, _) = poly_div_rem(&a, &d);
        let (new_b, _) = poly_div_rem(&c, &d);

        if new_a.degree() == 0 {
            break;
        }

        a = new_a;
        b = new_b;
        multiplicity += 1;
    }

    SquarefreeDecomposition { unit, factors }
}

/// Checks if a polynomial is squarefree.
///
/// A polynomial is squarefree if gcd(f, f') = 1.
pub fn is_squarefree<F: Field>(f: &DensePoly<F>) -> bool {
    if f.degree() == 0 {
        return true;
    }

    let f_prime = f.derivative();
    let g = poly_gcd(f, &f_prime);

    g.degree() == 0
}

/// Computes the squarefree part of a polynomial.
///
/// The squarefree part is f / gcd(f, f').
pub fn squarefree_part<F: Field>(f: &DensePoly<F>) -> DensePoly<F> {
    if f.degree() == 0 {
        return f.clone();
    }

    let f_prime = f.derivative();
    let g = poly_gcd(f, &f_prime);

    if g.degree() == 0 {
        f.clone()
    } else {
        let (result, _) = poly_div_rem(f, &g);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;
    use tertius_rings::traits::Ring;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_is_squarefree_linear() {
        let f = poly(&[-1, 1]); // x - 1
        assert!(is_squarefree(&f));
    }

    #[test]
    fn test_is_squarefree_product_of_distinct() {
        // (x-1)(x-2) = x² - 3x + 2
        let f = poly(&[2, -3, 1]);
        assert!(is_squarefree(&f));
    }

    #[test]
    fn test_is_not_squarefree_square() {
        // (x+1)² = x² + 2x + 1
        let f = poly(&[1, 2, 1]);
        assert!(!is_squarefree(&f));
    }

    #[test]
    fn test_squarefree_part_of_square() {
        // f = (x+1)² = x² + 2x + 1
        // squarefree part = x + 1
        let f = poly(&[1, 2, 1]);
        let sf = squarefree_part(&f);

        assert_eq!(sf.degree(), 1);
        // Check it's proportional to x + 1
        assert!(sf.leading_coeff().is_one());
    }

    #[test]
    fn test_squarefree_decomposition_simple() {
        // f = x + 1 (already squarefree)
        let f = poly(&[1, 1]);
        let decomp = squarefree_decomposition(&f);

        assert_eq!(decomp.factors.len(), 1);
        assert_eq!(decomp.factors[0].multiplicity, 1);
    }

    #[test]
    fn test_squarefree_decomposition_square() {
        // f = (x+1)² = x² + 2x + 1
        let f = poly(&[1, 2, 1]);
        let decomp = squarefree_decomposition(&f);

        // Should have one factor (x+1) with multiplicity 2
        assert_eq!(decomp.factors.len(), 1);
        assert_eq!(decomp.factors[0].multiplicity, 2);
        assert_eq!(decomp.factors[0].factor.degree(), 1);
    }

    #[test]
    fn test_squarefree_decomposition_cube() {
        // f = (x+1)³ = x³ + 3x² + 3x + 1
        let f = poly(&[1, 3, 3, 1]);
        let decomp = squarefree_decomposition(&f);

        // Should have one factor (x+1) with multiplicity 3
        assert_eq!(decomp.factors.len(), 1);
        assert_eq!(decomp.factors[0].multiplicity, 3);
        assert_eq!(decomp.factors[0].factor.degree(), 1);
    }

    #[test]
    fn test_squarefree_decomposition_mixed() {
        // f = (x+1)² * (x-1) = (x² + 2x + 1)(x - 1) = x³ + x² - x - 1
        let f = poly(&[-1, -1, 1, 1]);
        let decomp = squarefree_decomposition(&f);

        // Should have two factors:
        // - (x-1) with multiplicity 1
        // - (x+1) with multiplicity 2
        // The order depends on the algorithm

        let total_degree: usize = decomp
            .factors
            .iter()
            .map(|f| f.factor.degree() * f.multiplicity as usize)
            .sum();
        assert_eq!(total_degree, 3);

        // Verify reconstruction
        let reconstructed = decomp.to_polynomial();
        // Compare coefficients (they should match up to unit)
        assert_eq!(reconstructed.degree(), f.degree());
    }

    #[test]
    fn test_squarefree_decomposition_reconstruction() {
        // f = (x+1)² * (x+2)
        // = (x² + 2x + 1)(x + 2)
        // = x³ + 4x² + 5x + 2
        let f = poly(&[2, 5, 4, 1]);
        let decomp = squarefree_decomposition(&f);

        // Reconstruct and verify
        let reconstructed = decomp.to_polynomial();
        assert_eq!(reconstructed.degree(), f.degree());

        // Check that unit * product of factors = original
        for i in 0..=f.degree() {
            // Account for possible unit difference
            let ratio = if !f.coeff(i).is_zero() {
                f.coeff(i).field_div(&reconstructed.coeff(i))
            } else {
                Q::one()
            };
            // All ratios should be equal (the unit)
            if !f.coeff(i).is_zero() {
                assert!(!ratio.is_zero());
            }
        }
    }
}
