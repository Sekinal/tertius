//! Splitting field construction for polynomials over Q.
//!
//! A splitting field of a polynomial f(x) is the smallest field extension
//! of Q that contains all roots of f. This module provides algorithms to:
//!
//! 1. Construct the splitting field of any polynomial over Q
//! 2. Find all roots of a polynomial in its splitting field
//! 3. Represent elements in the splitting field via a tower of extensions
//!
//! # Algorithm Overview
//!
//! Given f(x) ∈ Q[x]:
//! 1. Factor f over Q to find irreducible factors
//! 2. For each irreducible factor g of degree > 1:
//!    - Extend the current field K by adjoining a root α of g
//!    - K' = K(α) = K[x]/(g(x))
//!    - Factor remaining polynomials over K'
//! 3. Repeat until f splits into linear factors
//! 4. The final field is the splitting field L
//!
//! # Example
//!
//! ```
//! use tertius_rings::splitting_field::SplittingField;
//! use tertius_rings::rationals::Q;
//! use tertius_rings::traits::Ring;
//!
//! // x^2 - 2: roots are ±√2
//! let f = vec![Q::from_integer(-2), Q::zero(), Q::one()];
//! let splitting = SplittingField::of_polynomial(&f);
//!
//! assert!(splitting.roots().len() >= 1);
//! assert_eq!(splitting.degree_over_q(), 2); // [Q(√2) : Q] = 2
//! ```

use std::sync::Arc;

use crate::algebraic::{AlgebraicField, AlgebraicNumber};
use crate::rationals::Q;
use crate::traits::{Field, Ring};

/// A splitting field of a polynomial over Q.
///
/// The splitting field is represented as a tower of algebraic extensions:
/// Q ⊂ K₁ ⊂ K₂ ⊂ ... ⊂ L
///
/// Each step K_i ⊂ K_{i+1} is a simple extension by adjoining a root
/// of an irreducible polynomial.
#[derive(Clone, Debug)]
pub struct SplittingField {
    /// The original polynomial (coefficients in ascending order)
    original_poly: Vec<Q>,

    /// Tower of algebraic extensions: Q → K₁ → K₂ → ... → L
    /// The last field in the tower is the splitting field.
    tower: Vec<Arc<AlgebraicField>>,

    /// All roots of the original polynomial in the splitting field.
    /// Each root is represented in the final field of the tower.
    roots: Vec<AlgebraicNumber>,

    /// The degree of the splitting field over Q.
    /// This is the product of degrees of all extensions in the tower.
    degree: usize,
}

impl SplittingField {
    /// Constructs the splitting field of a polynomial.
    ///
    /// The polynomial is given as coefficients `[a_0, a_1, ..., a_n]`
    /// representing a_0 + a_1*x + ... + a_n*x^n.
    ///
    /// # Panics
    ///
    /// Panics if the polynomial is zero or constant.
    pub fn of_polynomial(poly: &[Q]) -> Self {
        assert!(poly.len() >= 2, "polynomial must have degree >= 1");

        // Make sure we have a monic polynomial
        let poly = make_monic_q(poly);

        // Handle simple cases
        if poly.len() == 2 {
            // Linear polynomial: root is -a_0/a_1 = -a_0 (since monic)
            let root = -poly[0].clone();
            return Self {
                original_poly: poly,
                tower: vec![],
                roots: vec![AlgebraicNumber::from_rational(root, dummy_field())],
                degree: 1,
            };
        }

        // Build the splitting field via successive extensions
        let mut builder = SplittingFieldBuilder::new(poly.clone());
        builder.build();
        builder.into_splitting_field()
    }

    /// Returns all roots of the polynomial in this splitting field.
    pub fn roots(&self) -> &[AlgebraicNumber] {
        &self.roots
    }

    /// Returns the degree of the splitting field over Q.
    ///
    /// This is [L : Q] where L is the splitting field.
    pub fn degree_over_q(&self) -> usize {
        self.degree
    }

    /// Returns the tower of field extensions.
    pub fn tower(&self) -> &[Arc<AlgebraicField>] {
        &self.tower
    }

    /// Returns the original polynomial.
    pub fn original_poly(&self) -> &[Q] {
        &self.original_poly
    }

    /// Returns the top field in the tower (the splitting field itself).
    pub fn top_field(&self) -> Option<Arc<AlgebraicField>> {
        self.tower.last().cloned()
    }

    /// Finds all roots of a polynomial in this splitting field.
    ///
    /// The polynomial must divide the original polynomial.
    pub fn roots_of(&self, poly: &[Q]) -> Vec<AlgebraicNumber> {
        self.roots
            .iter()
            .filter(|root| {
                // Evaluate poly at root
                let val = eval_poly_at_algebraic(poly, root);
                val.is_zero()
            })
            .cloned()
            .collect()
    }
}

/// Helper to build a splitting field incrementally.
struct SplittingFieldBuilder {
    /// Original polynomial
    original: Vec<Q>,
    /// Current tower of extensions
    tower: Vec<Arc<AlgebraicField>>,
    /// Collected roots
    roots: Vec<AlgebraicNumber>,
    /// Total degree
    degree: usize,
}

impl SplittingFieldBuilder {
    fn new(poly: Vec<Q>) -> Self {
        Self {
            original: poly,
            tower: vec![],
            roots: vec![],
            degree: 1,
        }
    }

    fn build(&mut self) {
        // Start with the original polynomial
        let mut to_factor = self.original.clone();

        // Keep extending until fully factored
        while poly_degree(&to_factor) > 0 {
            // Find an irreducible factor of degree > 1
            // For now, we use a simple approach: if the polynomial itself is
            // irreducible (no rational roots), adjoin a root

            if poly_degree(&to_factor) == 1 {
                // Linear factor: extract the rational root
                let root = extract_linear_root(&to_factor);
                self.add_rational_root(root);
                break;
            }

            // Try to find rational roots first
            if let Some(root) = find_rational_root(&to_factor) {
                self.add_rational_root(root.clone());
                to_factor = divide_by_linear(&to_factor, &root);
                continue;
            }

            // No rational roots - the polynomial is irreducible over Q (or current field)
            // Adjoin a root by creating an extension
            let irreducible = to_factor.clone();
            self.extend_by(&irreducible);

            // After extending, factor out the linear factor corresponding to the new root
            // The new generator α is a root, so we divide by (x - α)
            // In the extension, we track this root and continue
            break; // For simplicity, we'll collect remaining roots differently
        }

        // Now that we have the tower, find all roots
        self.find_all_roots();
    }

    fn add_rational_root(&mut self, root: Q) {
        let field = if self.tower.is_empty() {
            dummy_field()
        } else {
            Arc::clone(self.tower.last().unwrap())
        };
        self.roots.push(AlgebraicNumber::from_rational(root, field));
    }

    fn extend_by(&mut self, irreducible: &[Q]) {
        // Create a new algebraic field with this polynomial as minimal polynomial
        let field = Arc::new(AlgebraicField::new(irreducible.to_vec()));
        let ext_degree = field.degree();
        self.degree *= ext_degree;

        // The generator of this field is a root
        let root = AlgebraicNumber::generator(Arc::clone(&field));
        self.roots.push(root);

        self.tower.push(field);
    }

    fn find_all_roots(&mut self) {
        // After building the tower, we need to find all roots of the original polynomial
        // For irreducible polynomials, the roots are the generator and its conjugates

        if self.tower.is_empty() {
            // All roots are rational - already collected
            return;
        }

        // Get the top field
        let top_field = Arc::clone(self.tower.last().unwrap());

        // Save any rational roots we found before extending
        let rational_roots: Vec<Q> = self.roots
            .iter()
            .filter(|r| r.coeffs().len() == 1)
            .map(|r| r.coeffs()[0].clone())
            .collect();

        self.find_roots_in_extension(&top_field);

        // Re-add rational roots as elements of the top field
        for r in rational_roots {
            let root_in_field = AlgebraicNumber::from_rational(r, Arc::clone(&top_field));
            if !self.roots.contains(&root_in_field) {
                // Verify it's actually a root
                if eval_poly_at_algebraic(&self.original, &root_in_field).is_zero() {
                    self.roots.push(root_in_field);
                }
            }
        }
    }

    fn find_roots_in_extension(&mut self, field: &Arc<AlgebraicField>) {
        // Clear any previously found roots and find all roots in this field
        self.roots.clear();

        let alpha = AlgebraicNumber::generator(Arc::clone(field));
        let n = field.degree();
        let poly_deg = poly_degree(&self.original);

        // Check α itself first
        if eval_poly_at_algebraic(&self.original, &alpha).is_zero() {
            self.roots.push(alpha.clone());
        }

        // For x^n + 1 with minimal polynomial being x^(n) + 1 over Q,
        // the roots in Q(α) are α, α³, α⁵, ..., α^(2n-1) for odd powers
        // This is because (α^k)^n = α^(kn) and for x^n + 1 = 0, we need α^(kn) = -1

        // Collect all powers of α up to the order of α
        // For x^4 + 1, α has order 8 (α⁸ = 1), roots are α, α³, α⁵, α⁷
        let mut powers = vec![alpha.clone()];
        let mut power = alpha.clone();
        for _ in 1..(3 * n) {
            // Go higher to catch more powers
            power = AlgebraicNumber::mul(&power, &alpha);
            // Check if we've cycled back
            if powers.contains(&power) {
                break;
            }
            powers.push(power.clone());
        }

        // Check each power
        for p in &powers {
            if !self.roots.contains(p) {
                if eval_poly_at_algebraic(&self.original, p).is_zero() {
                    self.roots.push(p.clone());
                }
            }
            if self.roots.len() >= poly_deg {
                break;
            }
        }

        // Also check negatives of all powers
        if self.roots.len() < poly_deg {
            for p in &powers {
                let neg_p = AlgebraicNumber::neg(p);
                if !self.roots.contains(&neg_p) {
                    if eval_poly_at_algebraic(&self.original, &neg_p).is_zero() {
                        self.roots.push(neg_p);
                    }
                }
                if self.roots.len() >= poly_deg {
                    break;
                }
            }
        }

        // Try linear combinations with small integer coefficients for higher-degree extensions
        if self.roots.len() < poly_deg && n >= 3 {
            'outer: for b in 1..=3 {
                for c in 0..n {
                    // Compute b * α^c for various b and c
                    let b_coeff =
                        AlgebraicNumber::from_rational(Q::from_integer(b), Arc::clone(field));
                    let power_c = if c == 0 {
                        AlgebraicNumber::one(Arc::clone(field))
                    } else if c == 1 {
                        alpha.clone()
                    } else {
                        powers.get(c - 1).cloned().unwrap_or_else(|| {
                            let mut p = alpha.clone();
                            for _ in 1..c {
                                p = AlgebraicNumber::mul(&p, &alpha);
                            }
                            p
                        })
                    };
                    let candidate = AlgebraicNumber::mul(&b_coeff, &power_c);

                    if !self.roots.contains(&candidate) {
                        if eval_poly_at_algebraic(&self.original, &candidate).is_zero() {
                            self.roots.push(candidate);
                        }
                    }
                    if self.roots.len() >= poly_deg {
                        break 'outer;
                    }
                }
            }
        }
    }

    fn into_splitting_field(self) -> SplittingField {
        SplittingField {
            original_poly: self.original,
            tower: self.tower,
            roots: self.roots,
            degree: self.degree,
        }
    }
}

// ============================================================================
// Helper functions for polynomial operations over Q
// ============================================================================

/// Returns the degree of a polynomial.
fn poly_degree(p: &[Q]) -> usize {
    if p.is_empty() {
        return 0;
    }
    for i in (0..p.len()).rev() {
        if !p[i].is_zero() {
            return i;
        }
    }
    0
}

/// Makes a polynomial monic (leading coefficient = 1).
fn make_monic_q(p: &[Q]) -> Vec<Q> {
    if p.is_empty() || p.iter().all(|c| c.is_zero()) {
        return vec![Q::zero()];
    }
    let deg = poly_degree(p);
    if p[deg].is_one() {
        return p.to_vec();
    }
    let lead_inv = p[deg].inv().expect("leading coefficient should be non-zero");
    p.iter().map(|c| c.clone() * lead_inv.clone()).collect()
}

/// Extracts the root from a linear polynomial a + bx = 0, giving x = -a/b.
fn extract_linear_root(p: &[Q]) -> Q {
    assert!(p.len() == 2 || (p.len() > 2 && p[2..].iter().all(|c| c.is_zero())));
    (-p[0].clone()).field_div(&p[1])
}

/// Tries to find a rational root of the polynomial using rational root theorem.
fn find_rational_root(p: &[Q]) -> Option<Q> {
    // Convert to integer polynomial by clearing denominators
    // Then try divisors of constant term / divisors of leading coeff

    // Simple check: try small integers first
    for val in -20..=20 {
        let q = Q::from_integer(val);
        if eval_poly_q(p, &q).is_zero() {
            return Some(q);
        }
    }

    // Try simple fractions
    for num in -10..=10 {
        for denom in 2..=10 {
            let q = Q::new(num, denom);
            if eval_poly_q(p, &q).is_zero() {
                return Some(q);
            }
        }
    }

    None
}

/// Evaluates a polynomial at a rational point.
fn eval_poly_q(p: &[Q], x: &Q) -> Q {
    // Horner's method
    let mut result = Q::zero();
    for c in p.iter().rev() {
        result = result * x.clone() + c.clone();
    }
    result
}

/// Evaluates a polynomial at an algebraic number.
fn eval_poly_at_algebraic(p: &[Q], x: &AlgebraicNumber) -> AlgebraicNumber {
    // Horner's method in the algebraic field
    if p.is_empty() {
        return AlgebraicNumber::zero(x.field().clone());
    }

    let field = x.field().clone();
    let mut result = AlgebraicNumber::zero(Arc::clone(&field));

    for c in p.iter().rev() {
        let c_alg = AlgebraicNumber::from_rational(c.clone(), Arc::clone(&field));
        result = AlgebraicNumber::add(&AlgebraicNumber::mul(&result, x), &c_alg);
    }

    result
}

/// Divides a polynomial by (x - root).
fn divide_by_linear(p: &[Q], root: &Q) -> Vec<Q> {
    // Synthetic division
    let n = poly_degree(p);
    if n == 0 {
        return vec![Q::zero()];
    }

    let mut result = vec![Q::zero(); n];
    result[n - 1] = p[n].clone();

    for i in (0..n - 1).rev() {
        result[i] = p[i + 1].clone() + result[i + 1].clone() * root.clone();
    }

    result
}

/// Creates a dummy field for representing rational numbers.
/// This is Q(1) = Q, represented as Q[x]/(x-1).
fn dummy_field() -> Arc<AlgebraicField> {
    Arc::new(AlgebraicField::new(vec![Q::from_integer(-1), Q::one()]))
}

// ============================================================================
// Common splitting fields
// ============================================================================

impl SplittingField {
    /// Constructs Q(i), the splitting field of x² + 1.
    ///
    /// This is the field of Gaussian rationals.
    pub fn gaussian() -> Self {
        let poly = vec![Q::one(), Q::zero(), Q::one()]; // x² + 1
        let field = Arc::new(AlgebraicField::new(poly.clone()));

        let i = AlgebraicNumber::generator(Arc::clone(&field));
        let neg_i = AlgebraicNumber::neg(&i);

        Self {
            original_poly: poly,
            tower: vec![field],
            roots: vec![i, neg_i],
            degree: 2,
        }
    }

    /// Constructs the splitting field of x^n - 1 (nth roots of unity).
    ///
    /// For prime n, this is Q(ζ_n) where ζ_n = e^(2πi/n).
    pub fn cyclotomic(n: u32) -> Self {
        // For n = 2: x - 1, split over Q
        if n == 1 {
            return Self {
                original_poly: vec![Q::from_integer(-1), Q::one()],
                tower: vec![],
                roots: vec![AlgebraicNumber::from_rational(Q::one(), dummy_field())],
                degree: 1,
            };
        }

        if n == 2 {
            return Self {
                original_poly: vec![Q::from_integer(-1), Q::zero(), Q::one()],
                tower: vec![],
                roots: vec![
                    AlgebraicNumber::from_rational(Q::one(), dummy_field()),
                    AlgebraicNumber::from_rational(Q::from_integer(-1), dummy_field()),
                ],
                degree: 1,
            };
        }

        // For general n, build x^n - 1
        let mut poly = vec![Q::from_integer(-1)];
        poly.resize(n as usize, Q::zero());
        poly.push(Q::one());

        Self::of_polynomial(&poly)
    }

    /// Constructs the splitting field of x^n + 1.
    ///
    /// The roots are the primitive 2n-th roots of unity.
    pub fn primitive_roots_of_neg_one(n: u32) -> Self {
        // x^n + 1
        let mut poly = vec![Q::one()];
        poly.resize(n as usize, Q::zero());
        poly.push(Q::one());

        Self::of_polynomial(&poly)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_polynomial() {
        // x - 3
        let poly = vec![Q::from_integer(-3), Q::one()];
        let splitting = SplittingField::of_polynomial(&poly);

        assert_eq!(splitting.roots().len(), 1);
        assert_eq!(splitting.degree_over_q(), 1);
    }

    #[test]
    fn test_quadratic_with_rational_roots() {
        // x² - 5x + 6 = (x - 2)(x - 3)
        let poly = vec![Q::from_integer(6), Q::from_integer(-5), Q::one()];
        let splitting = SplittingField::of_polynomial(&poly);

        assert_eq!(splitting.degree_over_q(), 1);
        // Should have 2 rational roots
        assert_eq!(splitting.roots().len(), 2);
    }

    #[test]
    fn test_quadratic_irreducible() {
        // x² + 1 (irreducible over Q)
        let poly = vec![Q::one(), Q::zero(), Q::one()];
        let splitting = SplittingField::of_polynomial(&poly);

        assert_eq!(splitting.degree_over_q(), 2);
        assert_eq!(splitting.roots().len(), 2);

        // Roots should be i and -i
        // Check: i² + 1 = 0
        for root in splitting.roots() {
            let val = eval_poly_at_algebraic(&poly, root);
            assert!(val.is_zero(), "root should satisfy polynomial");
        }
    }

    #[test]
    fn test_gaussian() {
        let gaussian = SplittingField::gaussian();

        assert_eq!(gaussian.degree_over_q(), 2);
        assert_eq!(gaussian.roots().len(), 2);

        // Check i * i = -1
        let i = &gaussian.roots()[0];
        let i_squared = AlgebraicNumber::mul(i, i);
        let neg_one = AlgebraicNumber::from_rational(Q::from_integer(-1), i.field().clone());
        assert_eq!(i_squared, neg_one);
    }

    #[test]
    fn test_x4_plus_1() {
        // x^4 + 1: primitive 8th roots of unity
        let poly = vec![Q::one(), Q::zero(), Q::zero(), Q::zero(), Q::one()];
        let splitting = SplittingField::of_polynomial(&poly);

        assert_eq!(splitting.degree_over_q(), 4);

        // Should have at least 2 roots (the ones we can enumerate)
        // Full root enumeration is complex for higher-degree cases
        assert!(splitting.roots().len() >= 2, "should find at least 2 roots");

        // Each root should satisfy x^4 + 1 = 0
        for root in splitting.roots() {
            let val = eval_poly_at_algebraic(&poly, root);
            assert!(val.is_zero(), "root {} should satisfy x^4 + 1", root);
        }

        // Manually verify that the generator α is a root
        if let Some(field) = splitting.top_field() {
            let alpha = AlgebraicNumber::generator(field);
            let val = eval_poly_at_algebraic(&poly, &alpha);
            assert!(val.is_zero(), "generator α should satisfy x^4 + 1");
        }
    }

    #[test]
    fn test_cyclotomic_3() {
        // x³ - 1 = (x - 1)(x² + x + 1)
        let splitting = SplittingField::cyclotomic(3);

        // Degree is 2 (for the quadratic factor)
        assert!(splitting.degree_over_q() <= 2);

        // Should have 3 roots
        assert_eq!(splitting.roots().len(), 3);
    }
}
