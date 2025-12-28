//! Asymptotic comparison classes for limit computation.
//!
//! When computing limits at infinity, we need to compare the growth rates
//! of different subexpressions. This module provides the comparison operations
//! needed by the Gruntz algorithm.

use std::cmp::Ordering;

/// Asymptotic comparison result between two expressions as x → ∞.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComparisonClass {
    /// f(x) / g(x) → 0 (f grows slower than g)
    LessThan,
    /// f(x) / g(x) → c ≠ 0, ∞ (f and g grow at the same rate)
    Comparable,
    /// f(x) / g(x) → ∞ (f grows faster than g)
    GreaterThan,
    /// Cannot determine the comparison
    Unknown,
}

impl ComparisonClass {
    /// Returns true if self represents "less than" (slower growth).
    pub fn is_less(&self) -> bool {
        matches!(self, ComparisonClass::LessThan)
    }

    /// Returns true if self represents "greater than" (faster growth).
    pub fn is_greater(&self) -> bool {
        matches!(self, ComparisonClass::GreaterThan)
    }

    /// Returns true if the expressions are comparable (same growth rate).
    pub fn is_comparable(&self) -> bool {
        matches!(self, ComparisonClass::Comparable)
    }

    /// Reverses the comparison (swaps the two expressions).
    pub fn reverse(self) -> Self {
        match self {
            ComparisonClass::LessThan => ComparisonClass::GreaterThan,
            ComparisonClass::GreaterThan => ComparisonClass::LessThan,
            other => other,
        }
    }
}

impl From<ComparisonClass> for Option<Ordering> {
    fn from(class: ComparisonClass) -> Self {
        match class {
            ComparisonClass::LessThan => Some(Ordering::Less),
            ComparisonClass::Comparable => Some(Ordering::Equal),
            ComparisonClass::GreaterThan => Some(Ordering::Greater),
            ComparisonClass::Unknown => None,
        }
    }
}

/// Growth rate classification for the Gruntz algorithm.
#[derive(Clone, Debug, PartialEq)]
pub enum GrowthRate {
    /// Constant (bounded) as x → ∞
    Constant,
    /// Polynomial growth: O(x^n)
    Polynomial(i32),
    /// Exponential growth: O(e^(c*x^n))
    Exponential { coeff: f64, power: i32 },
    /// Super-exponential: O(e^(e^x)) or faster
    SuperExponential,
    /// Logarithmic: O(log(x)^n)
    Logarithmic(i32),
    /// Unknown growth rate
    Unknown,
}

impl GrowthRate {
    /// Creates a polynomial growth rate.
    pub fn poly(n: i32) -> Self {
        if n == 0 {
            GrowthRate::Constant
        } else {
            GrowthRate::Polynomial(n)
        }
    }

    /// Creates an exponential growth rate.
    pub fn exp(power: i32) -> Self {
        GrowthRate::Exponential {
            coeff: 1.0,
            power,
        }
    }

    /// Compares two growth rates.
    pub fn compare(&self, other: &Self) -> ComparisonClass {
        match (self, other) {
            // Constants are comparable
            (GrowthRate::Constant, GrowthRate::Constant) => ComparisonClass::Comparable,

            // Constant < everything else non-constant
            (GrowthRate::Constant, _) => ComparisonClass::LessThan,
            (_, GrowthRate::Constant) => ComparisonClass::GreaterThan,

            // Logarithmic < Polynomial
            (GrowthRate::Logarithmic(_), GrowthRate::Polynomial(_)) => ComparisonClass::LessThan,
            (GrowthRate::Polynomial(_), GrowthRate::Logarithmic(_)) => ComparisonClass::GreaterThan,

            // Logarithmic comparisons
            (GrowthRate::Logarithmic(a), GrowthRate::Logarithmic(b)) => {
                match a.cmp(b) {
                    Ordering::Less => ComparisonClass::LessThan,
                    Ordering::Equal => ComparisonClass::Comparable,
                    Ordering::Greater => ComparisonClass::GreaterThan,
                }
            }

            // Polynomial comparisons
            (GrowthRate::Polynomial(a), GrowthRate::Polynomial(b)) => {
                match a.cmp(b) {
                    Ordering::Less => ComparisonClass::LessThan,
                    Ordering::Equal => ComparisonClass::Comparable,
                    Ordering::Greater => ComparisonClass::GreaterThan,
                }
            }

            // Polynomial < Exponential
            (GrowthRate::Polynomial(_), GrowthRate::Exponential { .. }) => {
                ComparisonClass::LessThan
            }
            (GrowthRate::Exponential { .. }, GrowthRate::Polynomial(_)) => {
                ComparisonClass::GreaterThan
            }

            // Logarithmic < Exponential
            (GrowthRate::Logarithmic(_), GrowthRate::Exponential { .. }) => {
                ComparisonClass::LessThan
            }
            (GrowthRate::Exponential { .. }, GrowthRate::Logarithmic(_)) => {
                ComparisonClass::GreaterThan
            }

            // Exponential comparisons (by power)
            (
                GrowthRate::Exponential { power: p1, .. },
                GrowthRate::Exponential { power: p2, .. },
            ) => {
                match p1.cmp(p2) {
                    Ordering::Less => ComparisonClass::LessThan,
                    Ordering::Equal => ComparisonClass::Comparable,
                    Ordering::Greater => ComparisonClass::GreaterThan,
                }
            }

            // SuperExponential > Exponential
            (GrowthRate::Exponential { .. }, GrowthRate::SuperExponential) => {
                ComparisonClass::LessThan
            }
            (GrowthRate::SuperExponential, GrowthRate::Exponential { .. }) => {
                ComparisonClass::GreaterThan
            }

            // SuperExponential > Polynomial
            (GrowthRate::Polynomial(_), GrowthRate::SuperExponential) => ComparisonClass::LessThan,
            (GrowthRate::SuperExponential, GrowthRate::Polynomial(_)) => {
                ComparisonClass::GreaterThan
            }

            // SuperExponential > Logarithmic
            (GrowthRate::Logarithmic(_), GrowthRate::SuperExponential) => {
                ComparisonClass::LessThan
            }
            (GrowthRate::SuperExponential, GrowthRate::Logarithmic(_)) => {
                ComparisonClass::GreaterThan
            }

            // SuperExponential comparable to itself
            (GrowthRate::SuperExponential, GrowthRate::SuperExponential) => {
                ComparisonClass::Comparable
            }

            // Unknown cases
            (GrowthRate::Unknown, _) | (_, GrowthRate::Unknown) => ComparisonClass::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comparison_class() {
        assert!(ComparisonClass::LessThan.is_less());
        assert!(ComparisonClass::GreaterThan.is_greater());
        assert!(ComparisonClass::Comparable.is_comparable());

        assert_eq!(
            ComparisonClass::LessThan.reverse(),
            ComparisonClass::GreaterThan
        );
    }

    #[test]
    fn test_growth_rate_ordering() {
        // Constant < Polynomial
        assert_eq!(
            GrowthRate::Constant.compare(&GrowthRate::Polynomial(1)),
            ComparisonClass::LessThan
        );

        // Polynomial < Polynomial (by degree)
        assert_eq!(
            GrowthRate::Polynomial(1).compare(&GrowthRate::Polynomial(2)),
            ComparisonClass::LessThan
        );

        // Polynomial < Exponential
        assert_eq!(
            GrowthRate::Polynomial(100).compare(&GrowthRate::Exponential {
                coeff: 1.0,
                power: 1
            }),
            ComparisonClass::LessThan
        );

        // Exponential < SuperExponential
        assert_eq!(
            GrowthRate::Exponential {
                coeff: 1.0,
                power: 1
            }
            .compare(&GrowthRate::SuperExponential),
            ComparisonClass::LessThan
        );
    }

    #[test]
    fn test_logarithmic() {
        // Log < Polynomial
        assert_eq!(
            GrowthRate::Logarithmic(1).compare(&GrowthRate::Polynomial(1)),
            ComparisonClass::LessThan
        );

        // Log(1) < Log(2)
        assert_eq!(
            GrowthRate::Logarithmic(1).compare(&GrowthRate::Logarithmic(2)),
            ComparisonClass::LessThan
        );
    }
}
