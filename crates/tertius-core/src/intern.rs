//! Hash-consing and interning utilities.
//!
//! This module provides the interning infrastructure that ensures
//! structural uniqueness of expressions in the arena.

use hashbrown::HashMap;
use std::hash::Hash;

/// A generic interning table.
///
/// This maps values to unique IDs, ensuring each unique value
/// is stored exactly once.
#[derive(Debug)]
pub struct InternTable<T> {
    /// Maps values to their IDs.
    map: HashMap<T, u32>,
    /// Stores values by ID for reverse lookup.
    values: Vec<T>,
}

impl<T: Clone + Eq + Hash> Default for InternTable<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Clone + Eq + Hash> InternTable<T> {
    /// Creates a new empty interning table.
    #[must_use]
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            values: Vec::new(),
        }
    }

    /// Creates a table with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            map: HashMap::with_capacity(capacity),
            values: Vec::with_capacity(capacity),
        }
    }

    /// Interns a value, returning its unique ID.
    ///
    /// If the value already exists, returns the existing ID.
    /// Otherwise, assigns a new ID and stores the value.
    pub fn intern(&mut self, value: T) -> u32 {
        if let Some(&id) = self.map.get(&value) {
            return id;
        }

        let id = self.values.len() as u32;
        self.map.insert(value.clone(), id);
        self.values.push(value);
        id
    }

    /// Gets a value by its ID.
    #[must_use]
    pub fn get(&self, id: u32) -> Option<&T> {
        self.values.get(id as usize)
    }

    /// Gets the ID of a value, if it exists.
    #[must_use]
    pub fn get_id(&self, value: &T) -> Option<u32> {
        self.map.get(value).copied()
    }

    /// Returns the number of interned values.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns true if no values have been interned.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Returns an iterator over all interned values.
    pub fn iter(&self) -> impl Iterator<Item = (u32, &T)> {
        self.values.iter().enumerate().map(|(i, v)| (i as u32, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intern_table() {
        let mut table = InternTable::new();

        let id1 = table.intern("hello".to_string());
        let id2 = table.intern("world".to_string());
        let id3 = table.intern("hello".to_string());

        assert_eq!(id1, 0);
        assert_eq!(id2, 1);
        assert_eq!(id1, id3); // Same value, same ID

        assert_eq!(table.get(id1), Some(&"hello".to_string()));
        assert_eq!(table.get(id2), Some(&"world".to_string()));
        assert_eq!(table.len(), 2);
    }
}
