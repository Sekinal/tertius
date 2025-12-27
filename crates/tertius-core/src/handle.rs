//! Type-safe expression handles.
//!
//! Handles are 32-bit indices into the arena, providing a lightweight
//! alternative to pointers with better cache performance.

use std::fmt;

/// A handle to an expression in the arena.
///
/// This is a lightweight 32-bit index that can be copied freely.
/// Two handles are equal if and only if they point to the same
/// (structurally identical) expression, thanks to hash-consing.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct ExprHandle(u32);

impl ExprHandle {
    /// Creates a new handle from an index.
    ///
    /// This is primarily for internal use by the arena.
    #[must_use]
    pub const fn new(index: u32) -> Self {
        Self(index)
    }

    /// Returns the raw index of this handle.
    #[must_use]
    pub const fn index(self) -> u32 {
        self.0
    }
}

impl fmt::Debug for ExprHandle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Expr({})", self.0)
    }
}

impl fmt::Display for ExprHandle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "#{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_handle_equality() {
        let h1 = ExprHandle::new(42);
        let h2 = ExprHandle::new(42);
        let h3 = ExprHandle::new(43);

        assert_eq!(h1, h2);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_handle_size() {
        // Ensure handles are only 4 bytes
        assert_eq!(std::mem::size_of::<ExprHandle>(), 4);
    }
}
