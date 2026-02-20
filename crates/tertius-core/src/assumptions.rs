//! Symbol assumptions for domain-safe symbolic transformations.

use hashbrown::{HashMap, HashSet};

/// Supported symbolic assumptions.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Assumption {
    /// Symbol is real-valued.
    Real,
    /// Symbol is strictly positive.
    Positive,
    /// Symbol is non-zero.
    NonZero,
    /// Symbol is an integer.
    Integer,
}

/// A mutable set of assumptions keyed by symbol name.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AssumptionSet {
    by_symbol: HashMap<String, HashSet<Assumption>>,
}

impl AssumptionSet {
    /// Creates an empty assumption set.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an assumption for a symbol.
    pub fn assume(&mut self, symbol: &str, assumption: Assumption) {
        self.by_symbol
            .entry(symbol.to_string())
            .or_default()
            .insert(assumption);
    }

    /// Returns true when the symbol has the assumption.
    #[must_use]
    pub fn has(&self, symbol: &str, assumption: Assumption) -> bool {
        self.by_symbol
            .get(symbol)
            .is_some_and(|set| set.contains(&assumption))
    }

    /// Returns true when all provided symbols are assumed real.
    #[must_use]
    pub fn all_real<'a>(&self, symbols: impl IntoIterator<Item = &'a str>) -> bool {
        symbols
            .into_iter()
            .all(|name| self.has(name, Assumption::Real) || self.has(name, Assumption::Positive))
    }

    /// Returns true if a symbol has any recorded assumptions.
    #[must_use]
    pub fn has_any(&self, symbol: &str) -> bool {
        self.by_symbol
            .get(symbol)
            .is_some_and(|set| !set.is_empty())
    }
}
