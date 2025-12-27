//! Expression node types.
//!
//! This module defines the core expression types stored in the arena.

use smallvec::SmallVec;

use crate::handle::ExprHandle;

/// Unique identifier for a symbol.
pub type SymbolId = u32;

/// Unique identifier for a function.
pub type FunctionId = u32;

/// An expression node stored in the arena.
///
/// This enum represents all possible expression types. Each variant is
/// designed to be cache-friendly, using `SmallVec` for inline storage
/// of small argument lists.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ExprNode {
    // === Atoms ===
    /// A 64-bit integer literal.
    ///
    /// For larger integers, use `BigInteger` (not yet implemented).
    Integer(i64),

    /// A rational number (numerator, denominator).
    ///
    /// Invariant: denominator > 0, gcd(num, den) == 1.
    Rational(i64, u64),

    /// A symbolic variable.
    Symbol(SymbolId),

    // === Compound Expressions ===
    /// Sum of expressions: a + b + c + ...
    ///
    /// Invariant: at least 2 arguments.
    Add(SmallVec<[ExprHandle; 4]>),

    /// Product of expressions: a * b * c * ...
    ///
    /// Invariant: at least 2 arguments.
    Mul(SmallVec<[ExprHandle; 4]>),

    /// Power expression: base^exp.
    Pow {
        /// The base of the power.
        base: ExprHandle,
        /// The exponent.
        exp: ExprHandle,
    },

    /// Negation: -expr.
    Neg(ExprHandle),

    /// Division: numerator / denominator.
    Div {
        /// The numerator.
        num: ExprHandle,
        /// The denominator.
        den: ExprHandle,
    },

    // === Functions ===
    /// A function application: f(arg1, arg2, ...).
    Function {
        /// The function identifier.
        id: FunctionId,
        /// The arguments.
        args: SmallVec<[ExprHandle; 2]>,
    },
}

impl ExprNode {
    /// Returns true if this node is an atom (no children).
    #[must_use]
    pub fn is_atom(&self) -> bool {
        matches!(
            self,
            ExprNode::Integer(_) | ExprNode::Rational(_, _) | ExprNode::Symbol(_)
        )
    }

    /// Returns true if this node is a numeric literal.
    #[must_use]
    pub fn is_number(&self) -> bool {
        matches!(self, ExprNode::Integer(_) | ExprNode::Rational(_, _))
    }

    /// Returns true if this is the integer zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        matches!(self, ExprNode::Integer(0))
    }

    /// Returns true if this is the integer one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        matches!(self, ExprNode::Integer(1))
    }

    /// Returns the children of this node.
    #[must_use]
    pub fn children(&self) -> SmallVec<[ExprHandle; 4]> {
        match self {
            ExprNode::Integer(_) | ExprNode::Rational(_, _) | ExprNode::Symbol(_) => {
                SmallVec::new()
            }
            ExprNode::Add(args) | ExprNode::Mul(args) => args.clone(),
            ExprNode::Pow { base, exp } => smallvec::smallvec![*base, *exp],
            ExprNode::Neg(arg) => smallvec::smallvec![*arg],
            ExprNode::Div { num, den } => smallvec::smallvec![*num, *den],
            ExprNode::Function { args, .. } => {
                args.iter().copied().collect()
            }
        }
    }
}

/// Standard function identifiers.
pub mod functions {
    use super::FunctionId;

    /// Sine function.
    pub const SIN: FunctionId = 0;
    /// Cosine function.
    pub const COS: FunctionId = 1;
    /// Tangent function.
    pub const TAN: FunctionId = 2;
    /// Natural exponential.
    pub const EXP: FunctionId = 3;
    /// Natural logarithm.
    pub const LN: FunctionId = 4;
    /// Logarithm base 10.
    pub const LOG10: FunctionId = 5;
    /// Square root.
    pub const SQRT: FunctionId = 6;
    /// Absolute value.
    pub const ABS: FunctionId = 7;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_atom() {
        assert!(ExprNode::Integer(42).is_atom());
        assert!(ExprNode::Symbol(0).is_atom());
        assert!(!ExprNode::Neg(ExprHandle::new(0)).is_atom());
    }

    #[test]
    fn test_is_zero_one() {
        assert!(ExprNode::Integer(0).is_zero());
        assert!(!ExprNode::Integer(1).is_zero());
        assert!(ExprNode::Integer(1).is_one());
        assert!(!ExprNode::Integer(0).is_one());
    }
}
