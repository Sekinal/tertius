//! The expression language for egg-based simplification.
//!
//! This module defines the language understood by the e-graph,
//! mapping Tertius expressions to egg nodes.

use egg::{define_language, Id, Symbol};

define_language! {
    /// The symbolic expression language for Tertius.
    pub enum TertiusLang {
        // Numeric literals
        Num(i64),
        // Variables
        Symbol(Symbol),

        // Basic arithmetic
        "+" = Add([Id; 2]),
        "-" = Sub([Id; 2]),
        "*" = Mul([Id; 2]),
        "/" = Div([Id; 2]),
        "neg" = Neg(Id),
        "^" = Pow([Id; 2]),

        // Trigonometric functions
        "sin" = Sin(Id),
        "cos" = Cos(Id),
        "tan" = Tan(Id),
        "asin" = Asin(Id),
        "acos" = Acos(Id),
        "atan" = Atan(Id),

        // Exponential and logarithmic
        "exp" = Exp(Id),
        "ln" = Ln(Id),
        "log10" = Log10(Id),
        "sqrt" = Sqrt(Id),

        // Other functions
        "abs" = Abs(Id),
        "factorial" = Factorial(Id),
    }
}

impl TertiusLang {
    /// Returns true if this node is a number.
    pub fn is_num(&self) -> bool {
        matches!(self, TertiusLang::Num(_))
    }

    /// Returns true if this node is a symbol.
    pub fn is_symbol(&self) -> bool {
        matches!(self, TertiusLang::Symbol(_))
    }

    /// Extracts the numeric value if this is a number.
    pub fn as_num(&self) -> Option<i64> {
        match self {
            TertiusLang::Num(n) => Some(*n),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egg::RecExpr;

    #[test]
    fn test_parse_expr() {
        let expr: RecExpr<TertiusLang> = "(+ 1 2)".parse().unwrap();
        assert_eq!(expr.as_ref().len(), 3);
    }

    #[test]
    fn test_parse_complex() {
        let expr: RecExpr<TertiusLang> = "(* (+ x 1) (- x 1))".parse().unwrap();
        assert!(!expr.as_ref().is_empty());
    }
}
