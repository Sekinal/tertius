//! Arena allocator for expression storage.
//!
//! This module provides a contiguous memory arena for storing expression nodes,
//! enabling cache-friendly traversal and constant-time deallocation.

use hashbrown::HashMap;
use smallvec::SmallVec;

use crate::expr::ExprNode;
use crate::handle::ExprHandle;

/// The main arena for storing expressions.
///
/// All expressions are stored contiguously in a `Vec`, with hash-consing
/// ensuring each unique expression is stored exactly once.
#[derive(Debug, Default)]
pub struct ExprArena {
    /// Storage for all expression nodes.
    nodes: Vec<ExprNode>,
    /// Interning table: maps node content to its handle.
    intern_map: HashMap<ExprNode, ExprHandle>,
    /// Symbol table: maps symbol names to their IDs.
    symbols: HashMap<String, u32>,
    /// Reverse symbol table for display.
    symbol_names: Vec<String>,
}

impl ExprArena {
    /// Creates a new empty arena.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates an arena with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(capacity),
            intern_map: HashMap::with_capacity(capacity),
            symbols: HashMap::new(),
            symbol_names: Vec::new(),
        }
    }

    /// Interns an expression node, returning its handle.
    ///
    /// If an identical node already exists, returns the existing handle.
    /// Otherwise, allocates a new node and returns its handle.
    pub fn intern(&mut self, node: ExprNode) -> ExprHandle {
        if let Some(&handle) = self.intern_map.get(&node) {
            return handle;
        }

        let index = self.nodes.len();
        assert!(index < u32::MAX as usize, "Arena capacity exceeded");

        let handle = ExprHandle::new(index as u32);
        self.nodes.push(node.clone());
        self.intern_map.insert(node, handle);
        handle
    }

    /// Gets the node at the given handle.
    ///
    /// # Panics
    ///
    /// Panics if the handle is invalid.
    #[must_use]
    pub fn get(&self, handle: ExprHandle) -> &ExprNode {
        &self.nodes[handle.index() as usize]
    }

    /// Interns a symbol, returning its unique ID.
    pub fn intern_symbol(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.symbols.get(name) {
            return id;
        }

        let id = self.symbol_names.len() as u32;
        self.symbols.insert(name.to_string(), id);
        self.symbol_names.push(name.to_string());
        id
    }

    /// Gets the name of a symbol by its ID.
    #[must_use]
    pub fn symbol_name(&self, id: u32) -> Option<&str> {
        self.symbol_names.get(id as usize).map(String::as_str)
    }

    /// Returns the number of nodes in the arena.
    #[must_use]
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Returns true if the arena is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    // === Convenience constructors ===

    /// Creates an integer expression.
    pub fn integer(&mut self, value: i64) -> ExprHandle {
        self.intern(ExprNode::Integer(value))
    }

    /// Creates a symbol expression.
    pub fn symbol(&mut self, name: &str) -> ExprHandle {
        let id = self.intern_symbol(name);
        self.intern(ExprNode::Symbol(id))
    }

    /// Creates an addition expression.
    pub fn add(&mut self, args: impl Into<SmallVec<[ExprHandle; 4]>>) -> ExprHandle {
        let args = args.into();
        if args.len() == 1 {
            return args[0];
        }
        self.intern(ExprNode::Add(args))
    }

    /// Creates a multiplication expression.
    pub fn mul(&mut self, args: impl Into<SmallVec<[ExprHandle; 4]>>) -> ExprHandle {
        let args = args.into();
        if args.len() == 1 {
            return args[0];
        }
        self.intern(ExprNode::Mul(args))
    }

    /// Creates a power expression.
    pub fn pow(&mut self, base: ExprHandle, exp: ExprHandle) -> ExprHandle {
        self.intern(ExprNode::Pow { base, exp })
    }

    /// Creates a negation expression.
    pub fn neg(&mut self, arg: ExprHandle) -> ExprHandle {
        self.intern(ExprNode::Neg(arg))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arena_basic() {
        let mut arena = ExprArena::new();

        let x = arena.symbol("x");
        let y = arena.symbol("y");

        // Same symbol returns same handle
        let x2 = arena.symbol("x");
        assert_eq!(x, x2);

        // Different symbols are different
        assert_ne!(x, y);
    }

    #[test]
    fn test_hash_consing() {
        let mut arena = ExprArena::new();

        let x = arena.symbol("x");
        let one = arena.integer(1);

        // Create (x + 1) twice
        let sum1 = arena.add(smallvec::smallvec![x, one]);
        let sum2 = arena.add(smallvec::smallvec![x, one]);

        // Should be the same handle due to hash-consing
        assert_eq!(sum1, sum2);

        // Arena should only have 3 nodes: x, 1, (x + 1)
        assert_eq!(arena.len(), 3);
    }
}
