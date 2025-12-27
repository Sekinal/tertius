//! M4GB-style reductor cache for efficient reductions.
//!
//! The reductor store caches tail-reduced polynomials indexed by their
//! leading monomials, enabling fast lookups during reduction.

use crate::labeled_poly::LabeledPoly;
use crate::monomial::PackedMonomial;
use parking_lot::RwLock;
use rustc_hash::FxHashMap;
use std::sync::Arc;

/// A cache of polynomials for reduction, indexed by leading monomial.
///
/// This implements the M4GB caching strategy where we store polynomials
/// that have been fully reduced, enabling efficient reuse.
pub struct ReductorStore<R> {
    /// Polynomials indexed by leading monomial.
    reducers: RwLock<FxHashMap<MonomialKey, Arc<LabeledPoly<R>>>>,
    /// Monomials that have been fully processed (no further reductions needed).
    fully_reduced: RwLock<FxHashMap<MonomialKey, ()>>,
}

/// A hashable key for monomials.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct MonomialKey {
    exponents: [u16; 16],
    total_degree: u32,
}

impl From<&PackedMonomial> for MonomialKey {
    fn from(m: &PackedMonomial) -> Self {
        let mut exponents = [0u16; 16];
        let exps = m.exponents();
        let n = exps.len().min(16);
        exponents[..n].copy_from_slice(&exps[..n]);
        Self {
            exponents,
            total_degree: m.total_degree(),
        }
    }
}

impl<R: Clone + PartialEq + Send + Sync> ReductorStore<R> {
    /// Creates a new empty reductor store.
    pub fn new() -> Self {
        Self {
            reducers: RwLock::new(FxHashMap::default()),
            fully_reduced: RwLock::new(FxHashMap::default()),
        }
    }

    /// Inserts a polynomial into the cache.
    pub fn insert(&self, poly: LabeledPoly<R>) {
        if let Some(lm) = poly.leading_monomial() {
            let key = MonomialKey::from(lm);
            let mut reducers = self.reducers.write();
            reducers.insert(key, Arc::new(poly));
        }
    }

    /// Looks up a polynomial by leading monomial.
    pub fn get(&self, lm: &PackedMonomial) -> Option<Arc<LabeledPoly<R>>> {
        let key = MonomialKey::from(lm);
        let reducers = self.reducers.read();
        reducers.get(&key).cloned()
    }

    /// Finds a reductor for a given monomial.
    ///
    /// Returns a polynomial whose leading monomial divides the given monomial.
    pub fn find_reductor(&self, m: &PackedMonomial) -> Option<Arc<LabeledPoly<R>>> {
        let reducers = self.reducers.read();

        // First try exact match
        let key = MonomialKey::from(m);
        if let Some(poly) = reducers.get(&key) {
            return Some(poly.clone());
        }

        // Otherwise, search for a divisor
        // This is O(n) in the worst case; could be improved with a trie structure
        for (_, poly) in reducers.iter() {
            if let Some(lm) = poly.leading_monomial() {
                if m.is_divisible_by(lm) {
                    return Some(poly.clone());
                }
            }
        }

        None
    }

    /// Marks a monomial as fully reduced (no further work needed).
    pub fn mark_fully_reduced(&self, m: &PackedMonomial) {
        let key = MonomialKey::from(m);
        let mut fully_reduced = self.fully_reduced.write();
        fully_reduced.insert(key, ());
    }

    /// Checks if a monomial is marked as fully reduced.
    pub fn is_fully_reduced(&self, m: &PackedMonomial) -> bool {
        let key = MonomialKey::from(m);
        let fully_reduced = self.fully_reduced.read();
        fully_reduced.contains_key(&key)
    }

    /// Returns the number of cached reducers.
    pub fn len(&self) -> usize {
        self.reducers.read().len()
    }

    /// Returns true if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.reducers.read().is_empty()
    }

    /// Clears the cache.
    pub fn clear(&self) {
        self.reducers.write().clear();
        self.fully_reduced.write().clear();
    }

    /// Returns all reducers as a vector.
    pub fn all_reducers(&self) -> Vec<Arc<LabeledPoly<R>>> {
        self.reducers.read().values().cloned().collect()
    }
}

impl<R: Clone + PartialEq + Send + Sync> Default for ReductorStore<R> {
    fn default() -> Self {
        Self::new()
    }
}

/// A more efficient reductor store using a monomial trie.
///
/// This provides O(log n) lookups for finding divisors.
pub struct TrieReductorStore<R> {
    /// The root of the trie.
    root: RwLock<TrieNode<R>>,
    /// Number of stored polynomials.
    count: RwLock<usize>,
}

struct TrieNode<R> {
    /// Polynomial stored at this node (if any).
    poly: Option<Arc<LabeledPoly<R>>>,
    /// Children indexed by (variable index, exponent).
    children: FxHashMap<(u8, u16), Box<TrieNode<R>>>,
}

impl<R> Default for TrieNode<R> {
    fn default() -> Self {
        Self {
            poly: None,
            children: FxHashMap::default(),
        }
    }
}

impl<R: Clone + PartialEq + Send + Sync> TrieReductorStore<R> {
    /// Creates a new empty trie reductor store.
    pub fn new() -> Self {
        Self {
            root: RwLock::new(TrieNode::default()),
            count: RwLock::new(0),
        }
    }

    /// Inserts a polynomial into the trie.
    pub fn insert(&self, poly: LabeledPoly<R>) {
        if let Some(lm) = poly.leading_monomial() {
            let mut root = self.root.write();
            let mut node = &mut *root;

            // Navigate to the appropriate node
            for (i, &exp) in lm.exponents().iter().enumerate() {
                if exp > 0 {
                    let key = (i as u8, exp);
                    node = node
                        .children
                        .entry(key)
                        .or_insert_with(|| Box::new(TrieNode::default()));
                }
            }

            node.poly = Some(Arc::new(poly));
            *self.count.write() += 1;
        }
    }

    /// Finds a reductor whose leading monomial divides the given monomial.
    pub fn find_reductor(&self, m: &PackedMonomial) -> Option<Arc<LabeledPoly<R>>> {
        let root = self.root.read();
        self.find_reductor_recursive(&root, m, 0)
    }

    fn find_reductor_recursive(
        &self,
        node: &TrieNode<R>,
        m: &PackedMonomial,
        var_start: usize,
    ) -> Option<Arc<LabeledPoly<R>>> {
        // Check if this node has a polynomial
        if let Some(ref poly) = node.poly {
            if let Some(lm) = poly.leading_monomial() {
                if m.is_divisible_by(lm) {
                    return Some(poly.clone());
                }
            }
        }

        // Try children
        let exps = m.exponents();
        for var in var_start..exps.len() {
            let target_exp = exps[var];
            if target_exp == 0 {
                continue;
            }

            // Try all exponents <= target_exp for this variable
            for exp in 1..=target_exp {
                let key = (var as u8, exp);
                if let Some(child) = node.children.get(&key) {
                    if let Some(poly) = self.find_reductor_recursive(child, m, var + 1) {
                        return Some(poly);
                    }
                }
            }
        }

        None
    }

    /// Returns the number of stored polynomials.
    pub fn len(&self) -> usize {
        *self.count.read()
    }

    /// Returns true if empty.
    pub fn is_empty(&self) -> bool {
        *self.count.read() == 0
    }
}

impl<R: Clone + PartialEq + Send + Sync> Default for TrieReductorStore<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signature::Signature;

    fn make_poly(lm: &[u16]) -> LabeledPoly<i64> {
        let mono = PackedMonomial::new(lm);
        let terms = vec![(1i64, mono.clone())];
        LabeledPoly::new(terms, Signature::generator(0, lm.len()), mono.total_degree())
    }

    #[test]
    fn test_reductor_store_basic() {
        let store: ReductorStore<i64> = ReductorStore::new();

        let poly1 = make_poly(&[2, 0]); // x^2
        let poly2 = make_poly(&[0, 2]); // y^2

        store.insert(poly1);
        store.insert(poly2);

        assert_eq!(store.len(), 2);

        // Should find x^2 as a reductor for x^3
        let m = PackedMonomial::new(&[3, 0]);
        assert!(store.find_reductor(&m).is_some());

        // Should not find a reductor for xy
        let m = PackedMonomial::new(&[1, 1]);
        assert!(store.find_reductor(&m).is_none());
    }

    #[test]
    fn test_fully_reduced() {
        let store: ReductorStore<i64> = ReductorStore::new();
        let m = PackedMonomial::new(&[1, 2, 3]);

        assert!(!store.is_fully_reduced(&m));
        store.mark_fully_reduced(&m);
        assert!(store.is_fully_reduced(&m));
    }
}
