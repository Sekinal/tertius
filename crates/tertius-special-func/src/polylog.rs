//! Polylogarithms Li_n(x).
//!
//! The polylogarithm of order n is defined as:
//!
//! Li_n(x) = Σ_{k=1}^∞ x^k / k^n
//!
//! Special cases:
//! - Li_1(x) = -ln(1-x)
//! - Li_2(x) = dilogarithm
//! - Li_3(x) = trilogarithm
//!
//! # Key Properties
//!
//! - Li_n(1) = ζ(n) (Riemann zeta function)
//! - Li_n(0) = 0
//! - d/dx Li_n(x) = Li_{n-1}(x) / x
//!
//! # Integration
//!
//! Polylogarithms arise from integrals like:
//! - ∫ ln(x)/(1-x) dx = -Li_2(x)
//! - ∫ ln(1-x)/x dx = -Li_2(x)

use tertius_rings::rationals::Q;
use tertius_rings::traits::{Ring, Field};

/// A polylogarithm function Li_n(argument).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polylogarithm {
    /// The order n of the polylogarithm.
    pub order: i32,
    /// Symbolic representation of the argument (as a string for now).
    pub argument: String,
}

impl Polylogarithm {
    /// Creates a new polylogarithm Li_n(arg).
    pub fn new(order: i32, argument: &str) -> Self {
        Self {
            order,
            argument: argument.to_string(),
        }
    }

    /// Creates the dilogarithm Li_2(arg).
    pub fn dilog(argument: &str) -> Self {
        Self::new(2, argument)
    }

    /// Creates the trilogarithm Li_3(arg).
    pub fn trilog(argument: &str) -> Self {
        Self::new(3, argument)
    }

    /// Returns the derivative d/dx Li_n(x) = Li_{n-1}(x) / x.
    pub fn derivative(&self) -> String {
        if self.order <= 1 {
            format!("-1/((1-{}) * {})", self.argument, self.argument)
        } else {
            format!("Li_{}({}) / {}", self.order - 1, self.argument, self.argument)
        }
    }

    /// Returns the functional equation Li_n(x) + Li_n(1-x) relation.
    pub fn functional_equation(&self) -> Option<String> {
        if self.order == 2 {
            // Li_2(x) + Li_2(1-x) = π²/6 - ln(x)ln(1-x)
            Some(format!(
                "Li_2({0}) + Li_2(1-{0}) = π²/6 - ln({0})*ln(1-{0})",
                self.argument
            ))
        } else {
            None
        }
    }
}

impl std::fmt::Display for Polylogarithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.order {
            2 => write!(f, "Li₂({})", self.argument),
            3 => write!(f, "Li₃({})", self.argument),
            n => write!(f, "Li_{}({})", n, self.argument),
        }
    }
}

/// Creates a dilogarithm Li_2(arg).
pub fn dilog(arg: &str) -> Polylogarithm {
    Polylogarithm::dilog(arg)
}

/// Creates a trilogarithm Li_3(arg).
pub fn trilog(arg: &str) -> Polylogarithm {
    Polylogarithm::trilog(arg)
}

/// Evaluates the polylogarithm series to n terms for |x| < 1.
pub fn polylog_series(order: i32, x: Q, num_terms: usize) -> Q {
    let mut result = Q::from_integer(0);
    let mut power = x.clone();

    for k in 1..=num_terms {
        let k_q = Q::from_integer(k as i64);
        let k_pow_n = k_q.clone().pow(order as u32);

        // term = x^k / k^n
        let term = power.clone() * k_pow_n.inv().unwrap();
        result = result + term;

        power = power * x.clone();
    }

    result
}

/// Standard dilogarithm values.
pub struct DilogValues;

impl DilogValues {
    /// Li_2(0) = 0
    pub fn at_zero() -> Q {
        Q::from_integer(0)
    }

    /// Li_2(1) = π²/6 (returns None since π is irrational)
    pub fn at_one() -> Option<Q> {
        None // π²/6 is irrational
    }

    /// Li_2(1/2) = π²/12 - (ln 2)²/2 (returns None)
    pub fn at_half() -> Option<Q> {
        None // irrational
    }

    /// Li_2(-1) = -π²/12
    pub fn at_neg_one() -> Option<Q> {
        None // irrational
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polylog_creation() {
        let li2 = Polylogarithm::dilog("x");
        assert_eq!(li2.order, 2);
        assert_eq!(li2.argument, "x");
    }

    #[test]
    fn test_polylog_display() {
        let li2 = dilog("x");
        assert_eq!(format!("{}", li2), "Li₂(x)");

        let li3 = trilog("1-x");
        assert_eq!(format!("{}", li3), "Li₃(1-x)");
    }

    #[test]
    fn test_polylog_derivative() {
        let li3 = trilog("x");
        let deriv = li3.derivative();
        assert!(deriv.contains("Li_2"));
    }

    #[test]
    fn test_dilog_functional_equation() {
        let li2 = dilog("x");
        let eq = li2.functional_equation();
        assert!(eq.is_some());
        assert!(eq.unwrap().contains("π²/6"));
    }

    #[test]
    fn test_polylog_series() {
        // Li_2(1/2) ≈ 0.582...
        // Series: 1/2 + 1/8 + 1/27·(1/8) + ...
        let half = Q::new(1, 2);
        let approx = polylog_series(2, half, 10);

        // Check it's positive and less than 1
        assert!(approx > Q::from_integer(0));
        assert!(approx < Q::from_integer(1));
    }
}
