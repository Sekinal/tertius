//! FGLM algorithm for Gröbner basis conversion.
//!
//! The FGLM algorithm (Faugère, Gianni, Lazard, Mora) converts a Gröbner
//! basis from one term ordering to another. This is essential for solving
//! polynomial systems: compute in grevlex (fast), then convert to lex
//! (triangular form for solving).
//!
//! # Algorithm Overview
//!
//! Given a Gröbner basis G w.r.t. grevlex:
//! 1. Enumerate monomials in lex order
//! 2. Compute normal forms w.r.t. G
//! 3. Detect linear dependencies → new basis elements
//!
//! Complexity: O(nD³) where D is the dimension of K[x]/I

use crate::quotient::{QuotientElement, Staircase};
use rustc_hash::FxHashMap;
use std::ops::Neg;
use tertius_groebner::labeled_poly::LabeledPoly;
use tertius_groebner::monomial::PackedMonomial;
use tertius_rings::traits::Field;

/// Result of the FGLM algorithm.
#[derive(Clone, Debug)]
pub struct FGLMResult<R> {
    /// Gröbner basis in lex order.
    pub lex_basis: Vec<Vec<(R, PackedMonomial)>>,
    /// Standard monomials of the quotient ring.
    pub staircase: Staircase,
    /// Dimension of the quotient ring (= number of solutions counting multiplicity).
    pub dimension: usize,
}

/// Incremental Gaussian elimination for linear dependency detection.
///
/// Maintains a row echelon form and allows testing if new vectors
/// are in the span of existing vectors.
struct IncrementalEchelon<R: Field + Clone> {
    /// Rows in row echelon form. Each row has (pivot_column, coefficients).
    rows: Vec<(usize, Vec<R>)>,
    /// Number of columns (dimension of vectors).
    num_cols: usize,
}

impl<R: Field + Clone> IncrementalEchelon<R> {
    fn new() -> Self {
        Self {
            rows: Vec::new(),
            num_cols: 0,
        }
    }

    /// Extends the dimension (adds more columns).
    fn extend_columns(&mut self, new_size: usize) {
        if new_size > self.num_cols {
            for (_, row) in &mut self.rows {
                row.resize(new_size, R::zero());
            }
            self.num_cols = new_size;
        }
    }

    /// Tries to add a new vector.
    ///
    /// Returns Some(coefficients) if the vector is linearly dependent,
    /// where coefficients express the vector as a linear combination of
    /// the original input vectors (in order of insertion).
    ///
    /// Returns None if the vector is linearly independent (and adds it).
    fn try_add(&mut self, mut vec: Vec<R>) -> Option<Vec<R>> {
        // Extend to current width
        vec.resize(self.num_cols.max(vec.len()), R::zero());
        self.extend_columns(vec.len());

        // Track how this vector is expressed in terms of original inputs
        // dependency_coeffs[i] = coefficient of the i-th original input
        let num_original = self.rows.len();
        let mut dependency_coeffs = vec![R::zero(); num_original];

        // Reduce the vector using existing rows
        for (i, (pivot_col, pivot_row)) in self.rows.iter().enumerate() {
            if *pivot_col < vec.len() && !vec[*pivot_col].is_zero() {
                let scale = vec[*pivot_col].clone();
                // vec -= scale * pivot_row
                for (j, v) in pivot_row.iter().enumerate() {
                    if j < vec.len() {
                        vec[j] = vec[j].clone() - scale.clone() * v.clone();
                    }
                }
                // Track the dependency
                dependency_coeffs[i] = scale;
            }
        }

        // Find the pivot (first non-zero entry)
        let pivot_col = vec.iter().position(|x| !x.is_zero());

        match pivot_col {
            Some(col) => {
                // Independent: normalize and add to echelon form
                let pivot_val = vec[col].clone();
                let pivot_inv = pivot_val.inv().expect("non-zero pivot");
                for v in &mut vec {
                    *v = v.clone() * pivot_inv.clone();
                }
                self.rows.push((col, vec));
                None
            }
            None => {
                // Dependent: the vector is zero after reduction
                // dependency_coeffs tells us how to express the original vector
                Some(dependency_coeffs)
            }
        }
    }
}

/// Converts a Gröbner basis from grevlex to lex ordering.
///
/// # Arguments
/// * `grevlex_basis` - Gröbner basis computed with grevlex ordering
/// * `num_vars` - Number of variables
///
/// # Returns
/// The Gröbner basis in lex ordering along with structural information.
pub fn fglm_convert<R>(grevlex_basis: &[LabeledPoly<R>], num_vars: usize) -> FGLMResult<R>
where
    R: Field + Clone + Send + Sync + Neg<Output = R> + std::fmt::Debug,
{
    #[cfg(test)]
    {
        eprintln!(
            "fglm_convert: starting with {} basis elements, {} vars",
            grevlex_basis.len(),
            num_vars
        );
        for (i, p) in grevlex_basis.iter().enumerate() {
            eprintln!("  basis[{}]: {} terms", i, p.num_terms());
            if let Some(lm) = p.leading_monomial() {
                eprintln!("    leading monomial: {:?}", lm);
            }
            for (c, m) in p.terms() {
                eprintln!("    term: ({:?}, {:?})", c, m);
            }
        }
    }

    if grevlex_basis.is_empty() || num_vars == 0 {
        return FGLMResult {
            lex_basis: vec![],
            staircase: Staircase::new(vec![PackedMonomial::one(num_vars)]),
            dimension: 1,
        };
    }

    // Build the quotient ring structure
    let mut standard_monomials: Vec<PackedMonomial> = Vec::new();
    let mut lex_basis_polys: Vec<Vec<(R, PackedMonomial)>> = Vec::new();

    // Incremental Gaussian elimination for linear independence
    let mut echelon: IncrementalEchelon<R> = IncrementalEchelon::new();

    // Track which monomials we've processed
    let mut processed: FxHashMap<PackedMonomial, bool> = FxHashMap::default(); // true = standard, false = reducible

    // Start with 1
    let one = PackedMonomial::one(num_vars);
    standard_monomials.push(one.clone());
    echelon.extend_columns(1);
    echelon.try_add(vec![R::one()]); // Add the first basis vector
    processed.insert(one.clone(), true);

    // Generate monomials in lex order using a worklist
    let mut worklist: std::collections::BinaryHeap<std::cmp::Reverse<LexMonomial>> =
        std::collections::BinaryHeap::new();

    // Add x_i for each variable
    for i in 0..num_vars {
        let m = PackedMonomial::var(i, num_vars);
        worklist.push(std::cmp::Reverse(LexMonomial(m)));
    }

    let mut iterations = 0;
    let max_iterations = 50000;
    let max_dimension = 10000;

    while let Some(std::cmp::Reverse(LexMonomial(current))) = worklist.pop() {
        iterations += 1;
        if iterations > max_iterations {
            break;
        }

        // Skip if we've already processed this monomial
        if processed.contains_key(&current) {
            continue;
        }

        #[cfg(test)]
        if iterations <= 10 {
            eprintln!("fglm_convert: iteration {}, processing {:?}", iterations, current);
        }

        // Compute normal form of current monomial w.r.t. grevlex basis
        let nf = normal_form_monomial(&current, grevlex_basis, num_vars);

        // Express normal form in terms of standard monomials
        let (nf_vector, is_complete) = express_in_basis(&nf, &standard_monomials);

        // If the normal form contains monomials not in the standard set,
        // then it's definitely independent
        if !is_complete {
            // Linearly independent: add to standard monomials
            processed.insert(current.clone(), true);
            standard_monomials.push(current.clone());

            // Extend echelon columns
            echelon.extend_columns(standard_monomials.len());

            // Add a row for this monomial (with 1 in its own position)
            let mut row = vec![R::zero(); standard_monomials.len()];
            row[standard_monomials.len() - 1] = R::one();
            echelon.try_add(row);

            // Add successors to worklist (multiply by each variable)
            for i in 0..num_vars {
                let successor = current.mul(&PackedMonomial::var(i, num_vars));
                if !processed.contains_key(&successor) {
                    worklist.push(std::cmp::Reverse(LexMonomial(successor)));
                }
            }
            continue;
        }

        // Try to add to the echelon form
        match echelon.try_add(nf_vector.clone()) {
            Some(coeffs) => {
                // Found a dependency: current - sum(coeffs[i] * standard[i]) ∈ I
                // This gives us a new lex basis polynomial
                processed.insert(current.clone(), false);

                let mut poly_terms: Vec<(R, PackedMonomial)> = Vec::new();
                poly_terms.push((R::one(), current.clone()));

                for (i, c) in coeffs.iter().enumerate() {
                    if !c.is_zero() && i < standard_monomials.len() {
                        poly_terms.push((c.clone().neg(), standard_monomials[i].clone()));
                    }
                }

                // Sort by lex descending (leading term first)
                poly_terms.sort_by(|a, b| b.1.cmp_lex(&a.1));

                // Remove zero coefficients
                poly_terms.retain(|(c, _)| !c.is_zero());

                if !poly_terms.is_empty() {
                    lex_basis_polys.push(poly_terms);
                }
            }
            None => {
                // Linearly independent: add to standard monomials
                processed.insert(current.clone(), true);
                standard_monomials.push(current.clone());

                // Extend echelon columns
                echelon.extend_columns(standard_monomials.len());

                // Add successors to worklist (multiply by each variable)
                for i in 0..num_vars {
                    let successor = current.mul(&PackedMonomial::var(i, num_vars));
                    if !processed.contains_key(&successor) {
                        worklist.push(std::cmp::Reverse(LexMonomial(successor)));
                    }
                }
            }
        }

        // Bound on dimension for safety
        if standard_monomials.len() > max_dimension {
            break;
        }
    }

    let staircase = Staircase::new(standard_monomials);
    let dimension = staircase.dimension();

    FGLMResult {
        lex_basis: lex_basis_polys,
        staircase,
        dimension,
    }
}

/// Wrapper for lexicographic ordering in BinaryHeap.
#[derive(Clone, Debug, Eq, PartialEq)]
struct LexMonomial(PackedMonomial);

impl Ord for LexMonomial {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp_lex(&other.0)
    }
}

impl PartialOrd for LexMonomial {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Computes the normal form of a monomial w.r.t. a Gröbner basis.
fn normal_form_monomial<R>(
    m: &PackedMonomial,
    basis: &[LabeledPoly<R>],
    num_vars: usize,
) -> QuotientElement<R>
where
    R: Field + Clone + Neg<Output = R>,
{
    let mut result = QuotientElement::from_monomial(m.clone(), R::one());

    let mut max_iters = 10000;
    loop {
        max_iters -= 1;
        if max_iters == 0 {
            break;
        }

        // Find a monomial in result that is reducible
        let reducible = result
            .coeffs
            .iter()
            .find(|(mon, coeff)| {
                !coeff.is_zero()
                    && basis.iter().any(|p| {
                        p.leading_monomial()
                            .map(|lm| mon.is_divisible_by(lm))
                            .unwrap_or(false)
                    })
            })
            .map(|(m, c)| (m.clone(), c.clone()));

        let (red_mon, red_coeff) = match reducible {
            Some((m, c)) => (m, c),
            None => break, // Fully reduced
        };

        // Find the reducer
        let reducer = basis.iter().find(|p| {
            p.leading_monomial()
                .map(|lm| red_mon.is_divisible_by(lm))
                .unwrap_or(false)
        });

        let reducer = match reducer {
            Some(r) => r,
            None => break,
        };

        let lm = reducer.leading_monomial().unwrap();
        let lc = reducer.leading_coeff().unwrap();
        let quotient = red_mon.div(lm).unwrap();

        // result -= (red_coeff / lc) * quotient * reducer
        // But the leading term of the reducer cancels with red_mon, so we only
        // subtract the tail terms.
        let scale = red_coeff.clone() * lc.inv().expect("lc invertible");

        // Remove red_mon from result
        result.coeffs.remove(&red_mon);

        // Add (not subtract!) scaled tail terms of reducer * quotient
        // Since reducer = lc*lm + tail, and we want: result - scale*(lc*lm + tail)
        // = result - scale*lc*lm - scale*tail
        // But we already removed red_mon = scale*lc*lm (since scale = red_coeff/lc and red_mon = quotient*lm)
        // Wait, actually red_mon IS the leading monomial * quotient = lm * quotient = lm * (red_mon/lm) = red_mon
        // So result - (red_coeff/lc) * reducer:
        //   result - (red_coeff/lc) * (lc * lm + tail)
        // = result - red_coeff * lm * quotient/quotient - (red_coeff/lc) * tail
        // Hmm, let me think again...
        //
        // We have result containing term: red_coeff * red_mon
        // reducer = lc * lm + (tail terms)
        // red_mon = quotient * lm (so quotient = red_mon / lm)
        //
        // We want to replace red_coeff * red_mon with:
        //   red_coeff * red_mon - (red_coeff/lc) * quotient * reducer
        // = red_coeff * red_mon - (red_coeff/lc) * quotient * (lc * lm + tail)
        // = red_coeff * red_mon - red_coeff * quotient * lm - (red_coeff/lc) * quotient * tail
        // = red_coeff * red_mon - red_coeff * red_mon - (red_coeff/lc) * quotient * tail
        // = - (red_coeff/lc) * quotient * tail
        //
        // So we need to subtract (scale * quotient * each tail term)
        for (c, mon) in reducer.terms() {
            // Skip the leading term - it cancels with red_mon
            if *mon == *lm {
                continue;
            }
            let new_mon = mon.mul(&quotient);
            let contribution = scale.clone() * c.clone();
            result
                .coeffs
                .entry(new_mon)
                .and_modify(|v| *v = v.clone() - contribution.clone())
                .or_insert_with(|| contribution.neg());
        }

        // Clean up zeros
        result.coeffs.retain(|_, c| !c.is_zero());
    }

    result
}

/// Expresses a quotient element in terms of standard monomials.
///
/// Returns (vector, is_complete) where is_complete indicates whether
/// all monomials in the element are in the standard set.
fn express_in_basis<R: Field + Clone>(
    elem: &QuotientElement<R>,
    standard: &[PackedMonomial],
) -> (Vec<R>, bool) {
    let standard_set: std::collections::HashSet<_> = standard.iter().collect();

    // Check if all monomials in elem are in the standard set
    let is_complete = elem.coeffs.keys().all(|m| standard_set.contains(m));

    let mut result = vec![R::zero(); standard.len()];
    for (i, m) in standard.iter().enumerate() {
        result[i] = elem.coeff(m);
    }
    (result, is_complete)
}

/// Solves a zero-dimensional polynomial system.
///
/// This is the main entry point for system solving:
/// 1. Compute Gröbner basis in grevlex
/// 2. Convert to lex using FGLM
/// 3. Extract solutions from triangular form
pub fn solve_system<R>(
    generators: Vec<Vec<(R, PackedMonomial)>>,
    num_vars: usize,
) -> FGLMResult<R>
where
    R: Field + Clone + Send + Sync + Neg<Output = R> + std::fmt::Debug,
{
    use tertius_groebner::m5gb::groebner_basis;

    // Compute grevlex Gröbner basis
    #[cfg(test)]
    eprintln!("solve_system: computing Gröbner basis...");
    let grevlex_basis = groebner_basis(generators);
    #[cfg(test)]
    eprintln!(
        "solve_system: Gröbner basis has {} elements",
        grevlex_basis.len()
    );

    // Convert to lex
    fglm_convert(&grevlex_basis, num_vars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_groebner::m5gb::groebner_basis;
    use tertius_rings::finite_field::FiniteField;

    type GF101 = FiniteField<101>;

    fn ff(n: u64) -> GF101 {
        GF101::new(n % 101)
    }

    fn mono(exps: &[u16]) -> PackedMonomial {
        PackedMonomial::new(exps)
    }

    #[test]
    fn test_fglm_simple_linear() {
        // System: x + y - 2 = 0, x - y = 0 over GF(101)
        // Solution: x = 1, y = 1

        // x + y + 99 (since 99 = -2 mod 101)
        let f = vec![
            (ff(1), mono(&[1, 0])),
            (ff(1), mono(&[0, 1])),
            (ff(99), mono(&[0, 0])),
        ];

        // x - y = x + 100*y (since 100 = -1 mod 101)
        let g = vec![(ff(1), mono(&[1, 0])), (ff(100), mono(&[0, 1]))];

        // Compute GB directly and then fglm_convert
        let gb = groebner_basis(vec![f, g]);
        let result = fglm_convert(&gb, 2);

        // Should have dimension 1 (unique solution)
        assert_eq!(result.dimension, 1);
        // The lex basis should give us x - 1 and y - 1 (or equivalent)
        assert!(!result.lex_basis.is_empty());
    }

    #[test]
    fn test_fglm_quadratic() {
        // System: x^2 - 1 = 0, y - x = 0
        // Solutions: (1, 1) and (-1, -1) = (100, 100) in GF(101)

        // x^2 - 1 = x^2 + 100 (since 100 = -1 mod 101)
        let f = vec![(ff(1), mono(&[2, 0])), (ff(100), mono(&[0, 0]))];

        // y - x = y + 100*x
        let g = vec![(ff(1), mono(&[0, 1])), (ff(100), mono(&[1, 0]))];

        // First compute GB directly (like test_fglm_univariate does)
        let gb = groebner_basis(vec![f, g]);

        // Then call fglm_convert
        let result = fglm_convert(&gb, 2);

        // Should have dimension 2 (two solutions)
        assert_eq!(result.dimension, 2);
        assert!(!result.lex_basis.is_empty());
    }

    #[test]
    fn test_incremental_echelon_independent() {
        let mut echelon: IncrementalEchelon<GF101> = IncrementalEchelon::new();
        echelon.extend_columns(3);

        // Add [1, 0, 0]
        let result = echelon.try_add(vec![ff(1), ff(0), ff(0)]);
        assert!(result.is_none()); // Independent

        // Add [0, 1, 0]
        let result = echelon.try_add(vec![ff(0), ff(1), ff(0)]);
        assert!(result.is_none()); // Independent

        // Add [0, 0, 1]
        let result = echelon.try_add(vec![ff(0), ff(0), ff(1)]);
        assert!(result.is_none()); // Independent
    }

    #[test]
    fn test_incremental_echelon_dependent() {
        let mut echelon: IncrementalEchelon<GF101> = IncrementalEchelon::new();
        echelon.extend_columns(3);

        // Add [1, 0, 0]
        echelon.try_add(vec![ff(1), ff(0), ff(0)]);

        // Add [0, 1, 0]
        echelon.try_add(vec![ff(0), ff(1), ff(0)]);

        // Add [1, 1, 0] = [1,0,0] + [0,1,0] - should be dependent
        let result = echelon.try_add(vec![ff(1), ff(1), ff(0)]);
        assert!(result.is_some());

        let coeffs = result.unwrap();
        // Should express [1,1,0] as 1*[1,0,0] + 1*[0,1,0]
        assert_eq!(coeffs[0], ff(1));
        assert_eq!(coeffs[1], ff(1));
    }

    #[test]
    fn test_compare_lex() {
        let a = mono(&[2, 0]);
        let b = mono(&[1, 1]);

        // In lex: x^2 > xy (compare first component first)
        assert_eq!(a.cmp_lex(&b), std::cmp::Ordering::Greater);
    }

    #[test]
    fn test_fglm_with_complete_basis() {
        // Use a known-complete Gröbner basis directly
        // System: x^2 - 1 = 0, x - y = 0
        // Complete GB (in grevlex): {x^2 - 1, x - y, y^2 - 1}
        use tertius_groebner::labeled_poly::LabeledPoly;

        let num_vars = 2;

        // x^2 - 1 (sorted by grevlex: x^2 > 1)
        let p1 = LabeledPoly::from_generator(
            vec![(ff(1), mono(&[2, 0])), (ff(100), mono(&[0, 0]))],
            0,
            num_vars,
        );
        // x - y (in grevlex: x > y since they have same degree, compare last var first)
        let p2 = LabeledPoly::from_generator(
            vec![(ff(1), mono(&[1, 0])), (ff(100), mono(&[0, 1]))],
            1,
            num_vars,
        );
        // y^2 - 1
        let p3 = LabeledPoly::from_generator(
            vec![(ff(1), mono(&[0, 2])), (ff(100), mono(&[0, 0]))],
            2,
            num_vars,
        );

        let complete_gb = vec![p1, p2, p3];
        let result = fglm_convert(&complete_gb, 2);

        // Should have dimension 2 (two solutions: (1,1) and (-1,-1))
        assert_eq!(result.dimension, 2);
        assert!(!result.lex_basis.is_empty());
    }

    #[test]
    fn test_fglm_univariate() {
        // Simple univariate: x^2 - 1 = 0 in one variable
        // Standard monomials: {1, x}

        let f = vec![(ff(1), mono(&[2])), (ff(100), mono(&[0]))]; // x^2 - 1

        let gb = groebner_basis(vec![f]);
        let result = fglm_convert(&gb, 1);

        // Dimension should be 2 (two roots: 1 and -1)
        assert_eq!(result.dimension, 2);
        assert_eq!(result.staircase.monomials.len(), 2);
    }

    #[test]
    fn test_fglm_cyclic2() {
        // cyclic-2: x + y = 0, x*y - 1 = 0
        // Solutions: x = 1, y = -1 and x = -1, y = 1

        // x + y
        let f1 = vec![(ff(1), mono(&[1, 0])), (ff(1), mono(&[0, 1]))];

        // xy - 1 = xy + 100
        let f2 = vec![(ff(1), mono(&[1, 1])), (ff(100), mono(&[0, 0]))];

        // Use direct groebner_basis + fglm_convert pattern
        let gb = groebner_basis(vec![f1, f2]);
        let result = fglm_convert(&gb, 2);

        // Should have dimension 2
        assert_eq!(result.dimension, 2);
        assert!(!result.lex_basis.is_empty());
    }

    /// Test just the groebner_basis call to isolate hanging issue
    #[test]
    fn test_groebner_only_bivariate() {
        // x^2 - 1, y - x
        let f = vec![(ff(1), mono(&[2, 0])), (ff(100), mono(&[0, 0]))];
        let g = vec![(ff(1), mono(&[0, 1])), (ff(100), mono(&[1, 0]))];

        eprintln!("About to call groebner_basis...");
        let gb = groebner_basis(vec![f, g]);
        eprintln!("groebner_basis returned {} elements", gb.len());

        assert!(gb.len() >= 3);
    }
}
