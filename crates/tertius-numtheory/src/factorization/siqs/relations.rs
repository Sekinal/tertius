//! Relation handling for SIQS.
//!
//! This module handles storing and processing relations for the linear algebra step.

use tertius_integers::Integer;

use super::Relation;

/// A relation database for SIQS.
pub struct RelationDatabase {
    /// Full relations (completely smooth).
    pub full: Vec<Relation>,
    /// Partial relations (smooth except for one large prime).
    pub partial: Vec<(u64, Relation)>,
    /// Factor base size.
    pub fb_size: usize,
}

impl RelationDatabase {
    /// Create a new relation database.
    pub fn new(fb_size: usize) -> Self {
        Self {
            full: Vec::new(),
            partial: Vec::new(),
            fb_size,
        }
    }

    /// Add a full relation.
    pub fn add_full(&mut self, rel: Relation) {
        self.full.push(rel);
    }

    /// Add a partial relation with a large prime.
    pub fn add_partial(&mut self, large_prime: u64, rel: Relation) {
        // Check if we can combine with an existing partial
        for i in 0..self.partial.len() {
            if self.partial[i].0 == large_prime {
                // Combine the two partials into a full relation
                let (_, other) = self.partial.remove(i);
                let combined = combine_partials(&rel, &other, large_prime);
                self.full.push(combined);
                return;
            }
        }

        self.partial.push((large_prime, rel));
    }

    /// Check if we have enough relations.
    pub fn has_enough(&self) -> bool {
        self.full.len() > self.fb_size + 5
    }

    /// Get all full relations.
    pub fn get_relations(&self) -> &[Relation] {
        &self.full
    }
}

/// Combine two partial relations that share a large prime.
fn combine_partials(a: &Relation, b: &Relation, _large_prime: u64) -> Relation {
    // The combined x is the product of the two x values
    // The exponents are XORed (since we're in GF(2))
    let x = a.x.clone() * b.x.clone();

    let mut exponents = a.exponents.clone();
    for (i, &e) in b.exponents.iter().enumerate() {
        if i < exponents.len() {
            exponents[i] ^= e;
        } else {
            exponents.push(e);
        }
    }

    let mut factorization = a.factorization.clone();
    for &(idx, exp) in &b.factorization {
        let mut found = false;
        for (fidx, fexp) in factorization.iter_mut() {
            if *fidx == idx {
                *fexp += exp;
                found = true;
                break;
            }
        }
        if !found {
            factorization.push((idx, exp));
        }
    }

    Relation {
        x,
        exponents,
        factorization,
    }
}
