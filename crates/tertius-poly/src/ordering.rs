//! Monomial orderings for polynomial operations.
//!
//! The choice of monomial ordering affects Gröbner basis computation
//! and other polynomial algorithms.

use std::cmp::Ordering;

use crate::monomial::{cmp_grevlex, cmp_grlex, cmp_lex, PackedMonomial};

/// A monomial ordering.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub enum MonomialOrder {
    /// Lexicographic order.
    ///
    /// x > y > z means x^a y^b z^c > x^d y^e z^f iff
    /// the first nonzero difference (a-d, b-e, c-f) is positive.
    Lex,

    /// Graded lexicographic order.
    ///
    /// First compares total degree, then uses lex as tiebreaker.
    Grlex,

    /// Graded reverse lexicographic order.
    ///
    /// First compares total degree, then uses reverse lex (last variable first)
    /// with the comparison reversed.
    #[default]
    Grevlex,
}

impl MonomialOrder {
    /// Compares two monomials according to this ordering.
    #[must_use]
    pub fn compare(&self, a: &PackedMonomial, b: &PackedMonomial, num_vars: usize) -> Ordering {
        match self {
            MonomialOrder::Lex => cmp_lex(a, b, num_vars),
            MonomialOrder::Grlex => cmp_grlex(a, b, num_vars),
            MonomialOrder::Grevlex => cmp_grevlex(a, b, num_vars),
        }
    }

    /// Returns a short name for the ordering.
    #[must_use]
    pub const fn name(&self) -> &'static str {
        match self {
            MonomialOrder::Lex => "lex",
            MonomialOrder::Grlex => "grlex",
            MonomialOrder::Grevlex => "grevlex",
        }
    }
}

impl std::fmt::Display for MonomialOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lex_order() {
        let order = MonomialOrder::Lex;

        let x = PackedMonomial::var(0, 2);
        let y = PackedMonomial::var(1, 2);
        let y2 = y.mul(&y);

        // x > y in lex
        assert_eq!(order.compare(&x, &y, 2), Ordering::Greater);

        // x > y^2 in lex (first variable dominates)
        assert_eq!(order.compare(&x, &y2, 2), Ordering::Greater);
    }

    #[test]
    fn test_grevlex_order() {
        let order = MonomialOrder::Grevlex;

        let x2 = PackedMonomial::from_exponents(&[2, 0]);
        let xy = PackedMonomial::from_exponents(&[1, 1]);
        let y2 = PackedMonomial::from_exponents(&[0, 2]);
        let x = PackedMonomial::var(0, 2);

        // Same degree: x^2 > xy > y^2
        assert_eq!(order.compare(&x2, &xy, 2), Ordering::Greater);
        assert_eq!(order.compare(&xy, &y2, 2), Ordering::Greater);

        // Higher degree wins: xy > x
        assert_eq!(order.compare(&xy, &x, 2), Ordering::Greater);
    }
}
