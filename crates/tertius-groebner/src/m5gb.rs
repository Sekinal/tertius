//! M5GB: Hybrid F5/M4GB Gröbner basis algorithm.
//!
//! This algorithm combines:
//! - F5 signatures for detecting useless pairs early
//! - M4GB-style caching for efficient reductions
//! - Parallel matrix reduction for batch processing

use crate::criteria::{is_rewritable_excluding, product_criterion, sugar_selection};
use crate::labeled_poly::LabeledPoly;
use crate::macaulay::{MacaulayMatrix, MonomialKey};
use crate::monomial::PackedMonomial;
use crate::parallel_reduce::{reduce_matrix, ReducedRow};
use crate::reductor_store::ReductorStore;
use crate::signature::{SignedPair, Signature};
use num_traits::Zero;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::ops::Neg;
use tertius_rings::traits::{Field, Ring};

/// Configuration for the M5GB algorithm.
#[derive(Clone, Debug)]
pub struct M5GBConfig {
    /// Maximum degree to process (0 = no limit).
    pub max_degree: u32,
    /// Whether to use F5 criteria.
    pub use_f5_criteria: bool,
    /// Whether to use the sugar selection strategy.
    pub use_sugar: bool,
    /// Maximum number of pairs to process per iteration.
    pub batch_size: usize,
}

impl Default for M5GBConfig {
    fn default() -> Self {
        Self {
            max_degree: 0,
            use_f5_criteria: true,
            use_sugar: true,
            batch_size: 100,
        }
    }
}

/// The M5GB Gröbner basis algorithm.
pub struct M5GB<R> {
    /// Current basis (labeled polynomials).
    basis: Vec<LabeledPoly<R>>,
    /// Pending pairs to process.
    pairs: Vec<SignedPair>,
    /// Reductor cache.
    cache: ReductorStore<R>,
    /// Processed pairs (for chain criterion).
    processed: Vec<(usize, usize)>,
    /// Configuration.
    config: M5GBConfig,
    /// Number of variables.
    num_vars: usize,
}

impl<R: Field + Clone + Send + Sync + Neg<Output = R>> M5GB<R> {
    /// Creates a new M5GB instance from generators.
    pub fn new(generators: Vec<Vec<(R, PackedMonomial)>>, config: M5GBConfig) -> Self {
        let num_vars = generators
            .first()
            .and_then(|g| g.first())
            .map(|(_, m)| m.num_vars())
            .unwrap_or(0);

        // Create labeled polynomials for each generator
        let mut basis: Vec<LabeledPoly<R>> = generators
            .into_iter()
            .enumerate()
            .filter_map(|(i, terms)| {
                if terms.is_empty() || terms.iter().all(|(c, _)| c.is_zero()) {
                    None
                } else {
                    // Sort terms by monomial descending
                    let mut terms = terms;
                    terms.sort_by(|a, b| b.1.cmp_grevlex(&a.1));
                    Some(LabeledPoly::from_generator(terms, i, num_vars))
                }
            })
            .collect();

        // Normalize: make all polynomials monic
        for poly in &mut basis {
            if let Some(lc) = poly.leading_coeff() {
                let inv = lc.inv().expect("leading coefficient must be invertible");
                let new_terms: Vec<_> = poly
                    .terms()
                    .iter()
                    .map(|(c, m)| (c.clone() * inv.clone(), m.clone()))
                    .collect();
                let sig = poly.signature.clone();
                let sugar = poly.sugar;
                *poly = LabeledPoly::new(new_terms, sig, sugar);
            }
        }

        // Generate initial pairs
        let mut pairs = Vec::new();
        for i in 0..basis.len() {
            for j in (i + 1)..basis.len() {
                if let (Some(lm_i), Some(lm_j)) =
                    (basis[i].leading_monomial(), basis[j].leading_monomial())
                {
                    // Skip if coprime (product criterion)
                    if product_criterion(lm_i, lm_j) {
                        continue;
                    }

                    let pair = SignedPair::new(
                        i,
                        j,
                        &basis[i].signature,
                        &basis[j].signature,
                        lm_i,
                        lm_j,
                        basis[i].sugar,
                        basis[j].sugar,
                    );
                    pairs.push(pair);
                }
            }
        }

        Self {
            basis,
            pairs,
            cache: ReductorStore::new(),
            processed: Vec::new(),
            config,
            num_vars,
        }
    }

    /// Computes the Gröbner basis.
    pub fn compute(mut self) -> Vec<LabeledPoly<R>> {
        while !self.pairs.is_empty() {
            self.step();
        }

        // Return the basis
        self.basis
    }

    /// Performs one step of the algorithm.
    fn step(&mut self) {
        // Select pairs to process
        let selected = if self.config.use_sugar {
            sugar_selection(&mut self.pairs)
        } else {
            // Take batch_size pairs
            let n = self.pairs.len().min(self.config.batch_size);
            self.pairs.drain(..n).collect()
        };

        if selected.is_empty() {
            return;
        }

        // Remove selected pairs from the list
        // (sugar_selection returns clones, so we need to remove them)
        if self.config.use_sugar {
            let min_sugar = selected.first().map(|p| p.sugar).unwrap_or(0);
            self.pairs.retain(|p| p.sugar != min_sugar);
        }

        // Filter pairs using F5 criteria
        // Note: We exclude the pair members themselves from the rewritability check
        let filtered: Vec<_> = if self.config.use_f5_criteria {
            selected
                .into_iter()
                .filter(|pair| {
                    !is_rewritable_excluding(&pair.signature, &self.basis, pair.i, pair.j)
                })
                .collect()
        } else {
            selected
        };

        if filtered.is_empty() {
            return;
        }

        // Compute S-polynomials and reduce
        let new_polys = self.reduce_spolys(&filtered);

        // Add new polynomials to basis
        for poly in new_polys {
            if !poly.is_zero() {
                self.add_to_basis(poly);
            }
        }

        // Mark pairs as processed
        for pair in filtered {
            let key = if pair.i < pair.j {
                (pair.i, pair.j)
            } else {
                (pair.j, pair.i)
            };
            self.processed.push(key);
        }
    }

    /// Reduces a batch of S-polynomials.
    fn reduce_spolys(&self, pairs: &[SignedPair]) -> Vec<LabeledPoly<R>> {
        // For each pair, compute the S-polynomial
        let spolys: Vec<_> = pairs
            .par_iter()
            .filter_map(|pair| self.compute_spoly(pair))
            .collect();

        // Reduce each S-polynomial by the basis
        spolys
            .into_par_iter()
            .map(|spoly| self.reduce_by_basis(spoly))
            .filter(|p| !p.is_zero())
            .collect()
    }

    /// Computes the S-polynomial for a pair.
    fn compute_spoly(&self, pair: &SignedPair) -> Option<LabeledPoly<R>> {
        let f = self.basis.get(pair.i)?;
        let g = self.basis.get(pair.j)?;

        let lm_f = f.leading_monomial()?;
        let lm_g = g.leading_monomial()?;

        let lcm = lm_f.lcm(lm_g);
        let mult_f = lcm.div(lm_f)?;
        let mult_g = lcm.div(lm_g)?;

        // S(f, g) = mult_f * f - mult_g * g
        let scaled_f = f.mul_monomial(&mult_f);
        let scaled_g = g.mul_monomial(&mult_g);

        // Subtract
        let result = self.subtract_polys(&scaled_f, &scaled_g, &pair.signature, pair.sugar);
        Some(result)
    }

    /// Subtracts two labeled polynomials.
    fn subtract_polys(
        &self,
        f: &LabeledPoly<R>,
        g: &LabeledPoly<R>,
        signature: &Signature,
        sugar: u32,
    ) -> LabeledPoly<R> {
        let mut terms: FxHashMap<MonomialKey, R> = FxHashMap::default();

        // Add terms from f
        for (c, m) in f.terms() {
            let key = MonomialKey::from(m);
            terms.insert(key, c.clone());
        }

        // Subtract terms from g
        for (c, m) in g.terms() {
            let key = MonomialKey::from(m);
            let neg_c = c.clone().neg();
            terms
                .entry(key)
                .and_modify(|v| *v = v.clone() + neg_c.clone())
                .or_insert(neg_c);
        }

        // Convert back to sorted list, filtering zeros
        let mut result: Vec<_> = terms
            .into_iter()
            .filter(|(_, c)| !c.is_zero())
            .map(|(key, c)| (c, PackedMonomial::new(&key.0)))
            .collect();

        result.sort_by(|a, b| b.1.cmp_grevlex(&a.1));

        LabeledPoly::new(result, signature.clone(), sugar)
    }

    /// Reduces a polynomial by the basis.
    fn reduce_by_basis(&self, mut poly: LabeledPoly<R>) -> LabeledPoly<R> {
        loop {
            let lm = match poly.leading_monomial() {
                Some(m) => m.clone(),
                None => break,
            };

            // Find a reducer
            let reducer = self.basis.iter().find(|g| {
                if let Some(lm_g) = g.leading_monomial() {
                    lm.is_divisible_by(lm_g)
                } else {
                    false
                }
            });

            let reducer = match reducer {
                Some(r) => r,
                None => break,
            };

            // Compute multiplier
            let lm_r = reducer.leading_monomial().unwrap();
            let mult = lm.div(lm_r).unwrap();

            // Compute the coefficient ratio
            let lc_poly = poly.leading_coeff().unwrap().clone();
            let lc_r = reducer.leading_coeff().unwrap().clone();
            let ratio = lc_poly * lc_r.inv().expect("leading coeff must be invertible");

            // Subtract ratio * mult * reducer from poly
            let scaled_reducer = reducer.mul_monomial(&mult);
            let mut terms: FxHashMap<MonomialKey, R> = FxHashMap::default();

            for (c, m) in poly.terms() {
                let key = MonomialKey::from(m);
                terms.insert(key, c.clone());
            }

            for (c, m) in scaled_reducer.terms() {
                let key = MonomialKey::from(m);
                let adjustment = (c.clone() * ratio.clone()).neg();
                terms
                    .entry(key)
                    .and_modify(|v| *v = v.clone() + adjustment.clone())
                    .or_insert(adjustment);
            }

            // Convert back
            let mut result: Vec<_> = terms
                .into_iter()
                .filter(|(_, c)| !c.is_zero())
                .map(|(key, c)| (c, PackedMonomial::new(&key.0)))
                .collect();

            result.sort_by(|a, b| b.1.cmp_grevlex(&a.1));

            let sig = poly.signature.clone();
            let sugar = poly.sugar;
            poly = LabeledPoly::new(result, sig, sugar);
        }

        // Make monic if non-zero
        if let Some(lc) = poly.leading_coeff() {
            let inv = lc.inv().expect("leading coefficient must be invertible");
            let new_terms: Vec<_> = poly
                .terms()
                .iter()
                .map(|(c, m)| (c.clone() * inv.clone(), m.clone()))
                .collect();
            let sig = poly.signature.clone();
            let sugar = poly.sugar;
            poly = LabeledPoly::new(new_terms, sig, sugar);
        }

        poly
    }

    /// Adds a polynomial to the basis and generates new pairs.
    fn add_to_basis(&mut self, poly: LabeledPoly<R>) {
        let new_idx = self.basis.len();
        let lm_new = match poly.leading_monomial() {
            Some(m) => m.clone(),
            None => return,
        };

        // Generate pairs with existing basis elements
        for (i, existing) in self.basis.iter().enumerate() {
            if let Some(lm_i) = existing.leading_monomial() {
                // Skip if coprime
                if product_criterion(lm_i, &lm_new) {
                    continue;
                }

                let pair = SignedPair::new(
                    i,
                    new_idx,
                    &existing.signature,
                    &poly.signature,
                    lm_i,
                    &lm_new,
                    existing.sugar,
                    poly.sugar,
                );

                // Check degree limit
                if self.config.max_degree > 0 && pair.lcm.total_degree() > self.config.max_degree {
                    continue;
                }

                self.pairs.push(pair);
            }
        }

        self.basis.push(poly);
    }

    /// Returns the current number of basis elements.
    pub fn basis_size(&self) -> usize {
        self.basis.len()
    }

    /// Returns the number of pending pairs.
    pub fn pairs_remaining(&self) -> usize {
        self.pairs.len()
    }
}

/// Convenience function to compute a Gröbner basis.
pub fn groebner_basis<R>(generators: Vec<Vec<(R, PackedMonomial)>>) -> Vec<LabeledPoly<R>>
where
    R: Field + Clone + Send + Sync + Neg<Output = R>,
{
    let m5gb = M5GB::new(generators, M5GBConfig::default());
    m5gb.compute()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::finite_field::FiniteField;

    type GF101 = FiniteField<101>;

    fn ff(n: u64) -> GF101 {
        GF101::new(n % 101)
    }

    fn mono(exps: &[u16]) -> PackedMonomial {
        PackedMonomial::new(exps)
    }

    #[test]
    fn test_m5gb_simple() {
        // f = x^2 - y
        let f = vec![(ff(1), mono(&[2, 0])), (ff(100), mono(&[0, 1]))]; // 100 = -1 mod 101

        // g = xy - 1
        let g = vec![(ff(1), mono(&[1, 1])), (ff(100), mono(&[0, 0]))];

        let gens = vec![f, g];
        let basis = groebner_basis(gens);

        // Should have a Gröbner basis
        assert!(!basis.is_empty());
    }

    #[test]
    fn test_m5gb_quadratic_with_linear() {
        // System that was failing in FGLM tests:
        // x^2 - 1 = 0, y - x = 0
        // Complete GB should be: {x^2 - 1, x - y, y^2 - 1} or equivalent

        // x^2 - 1
        let f = vec![(ff(1), mono(&[2, 0])), (ff(100), mono(&[0, 0]))];

        // y - x = y + 100*x
        let g = vec![(ff(1), mono(&[0, 1])), (ff(100), mono(&[1, 0]))];

        let gens = vec![f, g];
        let basis = groebner_basis(gens);

        // Should have 3 elements: x^2-1, x-y, y^2-1
        assert!(
            basis.len() >= 3,
            "Expected at least 3 basis elements, got {}",
            basis.len()
        );

        // Verify y^2 - 1 is in the basis (should have LM = y^2)
        let has_y2 = basis.iter().any(|p| {
            p.leading_monomial()
                .map(|m| m.exponent(0) == 0 && m.exponent(1) == 2)
                .unwrap_or(false)
        });
        assert!(has_y2, "Basis should contain y^2 - 1");
    }

    #[test]
    fn test_m5gb_linear() {
        // Simple linear system: x + y - 1, x - y - 1
        // Should give x - 1, y
        let f = vec![
            (ff(1), mono(&[1, 0])),
            (ff(1), mono(&[0, 1])),
            (ff(100), mono(&[0, 0])),
        ];
        let g = vec![
            (ff(1), mono(&[1, 0])),
            (ff(100), mono(&[0, 1])),
            (ff(100), mono(&[0, 0])),
        ];

        let gens = vec![f, g];
        let basis = groebner_basis(gens);

        // Basis should have 2 elements
        assert!(basis.len() >= 2);
    }

    #[test]
    fn test_already_groebner() {
        // Already a Gröbner basis: x, y
        let f = vec![(ff(1), mono(&[1, 0]))];
        let g = vec![(ff(1), mono(&[0, 1]))];

        let gens = vec![f, g];
        let basis = groebner_basis(gens);

        // Should still have 2 elements
        assert_eq!(basis.len(), 2);
    }
}
