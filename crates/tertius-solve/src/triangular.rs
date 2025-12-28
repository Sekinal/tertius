//! Triangular decomposition and solution extraction.
//!
//! A lex Gröbner basis for a zero-dimensional ideal has a triangular structure:
//! - g₁(x₁) - univariate in x₁
//! - g₂(x₁, x₂) - eliminates x₂
//! - ...
//! - gₙ(x₁, ..., xₙ) - all variables
//!
//! Solutions are found by solving g₁, then substituting into g₂, etc.

use num_traits::Zero;
use rustc_hash::FxHashMap;
use std::ops::Neg;
use tertius_groebner::monomial::PackedMonomial;
use tertius_rings::traits::Field;

/// A triangular system representation.
#[derive(Clone, Debug)]
pub struct TriangularSystem<R> {
    /// Polynomials in triangular form, indexed by elimination variable.
    /// polys[i] is the polynomial with leading variable x_i.
    pub polys: Vec<Option<Vec<(R, PackedMonomial)>>>,
    /// Number of variables.
    pub num_vars: usize,
}

impl<R: Field + Clone> TriangularSystem<R> {
    /// Extracts triangular structure from a lex Gröbner basis.
    ///
    /// For a zero-dimensional ideal, the lex basis has shape:
    /// - One polynomial with leading term x_n^{d_n}
    /// - One polynomial with leading term x_{n-1}^{d_{n-1}} (may involve x_n)
    /// - ...
    pub fn from_lex_basis(lex_basis: &[Vec<(R, PackedMonomial)>], num_vars: usize) -> Self {
        let mut polys: Vec<Option<Vec<(R, PackedMonomial)>>> = vec![None; num_vars];

        for poly in lex_basis {
            if poly.is_empty() {
                continue;
            }

            // Find the leading variable (first non-zero exponent in lex order)
            let leading_mon = &poly[0].1;
            let leading_var = find_leading_variable(leading_mon, num_vars);

            if let Some(var) = leading_var {
                // Keep the polynomial with the smallest degree in this variable
                let degree = leading_mon.exponent(var);
                let should_replace = polys[var].as_ref().map_or(true, |existing| {
                    let existing_degree = existing[0].1.exponent(var);
                    degree < existing_degree
                });

                if should_replace {
                    polys[var] = Some(poly.clone());
                }
            }
        }

        Self { polys, num_vars }
    }

    /// Checks if the system is complete (has a polynomial for each variable).
    pub fn is_complete(&self) -> bool {
        self.polys.iter().all(|p| p.is_some())
    }

    /// Returns the elimination polynomial for variable i.
    pub fn poly_for_var(&self, i: usize) -> Option<&Vec<(R, PackedMonomial)>> {
        self.polys.get(i).and_then(|p| p.as_ref())
    }
}

/// Finds the leading variable (first non-zero exponent in lex order).
fn find_leading_variable(m: &PackedMonomial, num_vars: usize) -> Option<usize> {
    for i in 0..num_vars {
        if m.exponent(i) > 0 {
            return Some(i);
        }
    }
    None
}

/// A solution point.
#[derive(Clone, Debug)]
pub struct Solution<R> {
    /// Values for each variable.
    pub values: Vec<R>,
}

impl<R: Clone> Solution<R> {
    /// Creates a new solution.
    pub fn new(values: Vec<R>) -> Self {
        Self { values }
    }

    /// Returns the value for variable i.
    pub fn get(&self, i: usize) -> Option<&R> {
        self.values.get(i)
    }
}

/// Solves a triangular system over a finite field.
///
/// For finite fields, we can enumerate all possible values and check.
/// This is practical for small fields.
pub fn solve_triangular<R>(system: &TriangularSystem<R>) -> Vec<Solution<R>>
where
    R: Field + Clone + Neg<Output = R> + FiniteFieldExt,
{
    let mut solutions = Vec::new();
    let num_vars = system.num_vars;

    if num_vars == 0 {
        return vec![Solution::new(vec![])];
    }

    // Find the first univariate polynomial (polynomial involving only one variable)
    // In a proper lex triangular system, we start from the first variable
    let first_var = find_univariate_var(system);

    if let Some(first_var) = first_var {
        if let Some(uni_poly) = system.poly_for_var(first_var) {
            // Find roots of the univariate polynomial
            let roots = find_roots_finite_field_univariate(uni_poly, first_var);

            for root in roots {
                // Try to extend this partial solution
                let mut partial = vec![R::zero(); num_vars];
                partial[first_var] = root;

                extend_solution_forward(system, &mut partial, first_var, &mut solutions);
            }
        }
    } else {
        // No univariate polynomial found - try brute force for first variable
        let mut partial = vec![R::zero(); num_vars];
        for val in R::all_elements() {
            partial[0] = val;
            extend_solution_forward(system, &mut partial, 0, &mut solutions);
        }
    }

    solutions
}

/// Finds the variable with a univariate polynomial.
fn find_univariate_var<R: Field + Clone>(system: &TriangularSystem<R>) -> Option<usize> {
    for (i, poly_opt) in system.polys.iter().enumerate() {
        if let Some(poly) = poly_opt {
            // Check if this polynomial only involves variable i
            if is_univariate(poly, i, system.num_vars) {
                return Some(i);
            }
        }
    }
    None
}

/// Checks if a polynomial is univariate in the given variable.
fn is_univariate<R>(poly: &[(R, PackedMonomial)], var: usize, num_vars: usize) -> bool {
    for (_, mono) in poly {
        for j in 0..num_vars {
            if j != var && mono.exponent(j) > 0 {
                return false;
            }
        }
    }
    true
}

/// Finds roots of a truly univariate polynomial.
fn find_roots_finite_field_univariate<R: Field + Clone + FiniteFieldExt>(
    poly: &[(R, PackedMonomial)],
    var: usize,
) -> Vec<R> {
    // Extract coefficients: coeff[i] is the coefficient of var^i
    let mut max_deg = 0u16;
    for (_, mono) in poly {
        max_deg = max_deg.max(mono.exponent(var));
    }

    let mut coeffs = vec![R::zero(); (max_deg + 1) as usize];
    for (c, mono) in poly {
        let exp = mono.exponent(var) as usize;
        coeffs[exp] = coeffs[exp].clone() + c.clone();
    }

    // Convert to (coeff, exp) format
    let uni: Vec<(R, u16)> = coeffs
        .into_iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, c)| (c, i as u16))
        .collect();

    find_roots_univariate(&uni)
}

/// Extends a solution forward (from smaller to larger variable indices).
fn extend_solution_forward<R>(
    system: &TriangularSystem<R>,
    partial: &mut Vec<R>,
    current_var: usize,
    solutions: &mut Vec<Solution<R>>,
) where
    R: Field + Clone + Neg<Output = R> + FiniteFieldExt,
{
    let num_vars = system.num_vars;

    if current_var >= num_vars - 1 {
        // All variables assigned - verify the solution
        if verify_solution(system, partial) {
            solutions.push(Solution::new(partial.clone()));
        }
        return;
    }

    let next_var = current_var + 1;

    if let Some(poly) = system.poly_for_var(next_var) {
        // Substitute known values and solve for next_var
        let substituted = substitute_partial(poly, partial, next_var);
        let roots = find_roots_univariate(&substituted);

        for root in roots {
            partial[next_var] = root;
            extend_solution_forward(system, partial, next_var, solutions);
        }
    } else {
        // No constraint on this variable - enumerate
        for val in R::all_elements() {
            partial[next_var] = val;
            extend_solution_forward(system, partial, next_var, solutions);
        }
    }
}

/// Substitutes known values into a polynomial.
fn substitute_partial<R: Field + Clone + Neg<Output = R>>(
    poly: &[(R, PackedMonomial)],
    values: &[R],
    target_var: usize,
) -> Vec<(R, u16)> {
    // Result: univariate polynomial in target_var
    let mut result: FxHashMap<u16, R> = FxHashMap::default();

    for (coeff, mono) in poly {
        let target_exp = mono.exponent(target_var);

        // Evaluate other variables at their known values
        let mut eval_coeff = coeff.clone();
        for (i, val) in values.iter().enumerate() {
            if i != target_var {
                let exp = mono.exponent(i);
                if exp > 0 {
                    eval_coeff = eval_coeff * val.pow(exp as u32);
                }
            }
        }

        result
            .entry(target_exp)
            .and_modify(|c| *c = c.clone() + eval_coeff.clone())
            .or_insert(eval_coeff);
    }

    // Convert to sorted vector (coeff, exp) pairs
    let mut terms: Vec<(R, u16)> = result
        .into_iter()
        .filter(|(_, c)| !c.is_zero())
        .map(|(exp, coeff)| (coeff, exp))
        .collect();
    terms.sort_by(|a, b| b.1.cmp(&a.1)); // Descending by degree
    terms
}

/// Finds roots of a univariate polynomial over a finite field.
fn find_roots_univariate<R: Field + Clone + FiniteFieldExt>(poly: &[(R, u16)]) -> Vec<R> {
    if poly.is_empty() {
        return R::all_elements(); // Zero polynomial - all elements are roots
    }

    let mut roots = Vec::new();

    for candidate in R::all_elements() {
        let value = evaluate_univariate(poly, &candidate);
        if value.is_zero() {
            roots.push(candidate);
        }
    }

    roots
}

/// Evaluates a univariate polynomial at a point.
fn evaluate_univariate<R: Field + Clone>(poly: &[(R, u16)], x: &R) -> R {
    let mut result = R::zero();
    for (coeff, exp) in poly {
        let term = coeff.clone() * x.pow(*exp as u32);
        result = result + term;
    }
    result
}

/// Verifies that a candidate solution satisfies all polynomials.
fn verify_solution<R: Field + Clone>(system: &TriangularSystem<R>, sol: &[R]) -> bool {
    for poly_opt in &system.polys {
        if let Some(poly) = poly_opt {
            let value = evaluate_poly(poly, sol);
            if !value.is_zero() {
                return false;
            }
        }
    }
    true
}

/// Evaluates a polynomial at a point.
fn evaluate_poly<R: Field + Clone>(poly: &[(R, PackedMonomial)], point: &[R]) -> R {
    let mut result = R::zero();
    for (coeff, mono) in poly {
        let mut term = coeff.clone();
        for (i, val) in point.iter().enumerate() {
            let exp = mono.exponent(i);
            if exp > 0 {
                term = term * val.pow(exp as u32);
            }
        }
        result = result + term;
    }
    result
}

/// Extension trait for finite field enumeration.
pub trait FiniteFieldExt: Sized {
    /// Returns an iterator over all field elements.
    fn all_elements() -> Vec<Self>;
}

// Implement for common finite field types
impl<const P: u64> FiniteFieldExt for tertius_rings::finite_field::FiniteField<P> {
    fn all_elements() -> Vec<Self> {
        (0..P).map(Self::new).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::finite_field::FiniteField;

    type GF7 = FiniteField<7>;

    fn ff(n: u64) -> GF7 {
        GF7::new(n % 7)
    }

    fn mono(exps: &[u16]) -> PackedMonomial {
        PackedMonomial::new(exps)
    }

    #[test]
    fn test_triangular_from_lex_basis() {
        // Simple lex basis: x - 1, y - 2
        let poly1 = vec![(ff(1), mono(&[1, 0])), (ff(6), mono(&[0, 0]))]; // x - 1 = x + 6
        let poly2 = vec![(ff(1), mono(&[0, 1])), (ff(5), mono(&[0, 0]))]; // y - 2 = y + 5

        let system = TriangularSystem::from_lex_basis(&[poly1, poly2], 2);

        assert!(system.is_complete());
        assert!(system.poly_for_var(0).is_some());
        assert!(system.poly_for_var(1).is_some());
    }

    #[test]
    fn test_solve_simple() {
        // System: x = 1, y = 2 in GF(7)
        let poly1 = vec![(ff(1), mono(&[1, 0])), (ff(6), mono(&[0, 0]))]; // x - 1
        let poly2 = vec![(ff(1), mono(&[0, 1])), (ff(5), mono(&[0, 0]))]; // y - 2

        let system = TriangularSystem::from_lex_basis(&[poly1, poly2], 2);
        let solutions = solve_triangular(&system);

        assert_eq!(solutions.len(), 1);
        assert_eq!(solutions[0].values, vec![ff(1), ff(2)]);
    }

    #[test]
    fn test_solve_quadratic() {
        // System: x^2 = 1, y = x in GF(7)
        // Solutions: (1, 1) and (6, 6) since 6 = -1 mod 7

        // x^2 - 1 = x^2 + 6
        let poly1 = vec![(ff(1), mono(&[2, 0])), (ff(6), mono(&[0, 0]))];
        // y - x
        let poly2 = vec![(ff(1), mono(&[0, 1])), (ff(6), mono(&[1, 0]))];

        let system = TriangularSystem::from_lex_basis(&[poly1, poly2], 2);
        let solutions = solve_triangular(&system);

        assert_eq!(solutions.len(), 2);

        // Check solutions are correct
        let sol_set: std::collections::HashSet<_> = solutions
            .iter()
            .map(|s| (s.values[0].value(), s.values[1].value()))
            .collect();

        assert!(sol_set.contains(&(1, 1)));
        assert!(sol_set.contains(&(6, 6)));
    }

    #[test]
    fn test_find_roots_univariate() {
        // x^2 - 1 in GF(7)
        let poly = vec![(ff(1), 2u16), (ff(6), 0u16)];
        let roots = find_roots_univariate(&poly);

        assert_eq!(roots.len(), 2);
        assert!(roots.contains(&ff(1)));
        assert!(roots.contains(&ff(6)));
    }
}
