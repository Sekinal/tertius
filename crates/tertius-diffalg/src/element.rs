//! Elements of a transcendental tower.
//!
//! An element of K(θ₁, ..., θₙ) can be represented as a rational function
//! in θₙ with coefficients in K(θ₁, ..., θₙ₋₁).
//!
//! This recursive structure mirrors the tower itself.

use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::traits::Field;

use crate::tower::TranscendentalTower;

/// An element of a transcendental tower.
///
/// At each level, an element is a rational function in the current
/// transcendental θ with coefficients from the level below.
#[derive(Clone, Debug)]
pub enum TowerElement<F: Field> {
    /// An element of the base field K = Q(x)
    Base(RationalFunction<F>),

    /// An element of an extension K(θ), represented as a polynomial
    /// in θ with coefficients that are TowerElements from the level below.
    Extension {
        /// The level of this element (1-indexed: 1 for K(θ₁), 2 for K(θ₁, θ₂), etc.)
        level: usize,
        /// Coefficients of the polynomial in θ (from constant term up)
        /// Each coefficient is an element of the previous level
        numerator_coeffs: Vec<TowerElement<F>>,
        /// Denominator coefficients (for rational function in θ)
        denominator_coeffs: Vec<TowerElement<F>>,
    },
}

impl<F: Field> TowerElement<F> {
    /// Creates a base field element from a rational function.
    pub fn from_base(rf: RationalFunction<F>) -> Self {
        TowerElement::Base(rf)
    }

    /// Creates a base field element from a polynomial.
    pub fn from_poly(p: DensePoly<F>) -> Self {
        TowerElement::Base(RationalFunction::from_poly(p))
    }

    /// Creates the zero element at the base level.
    pub fn zero() -> Self {
        TowerElement::Base(RationalFunction::zero())
    }

    /// Creates the one element at the base level.
    pub fn one() -> Self {
        TowerElement::Base(RationalFunction::one())
    }

    /// Returns the level of this element (0 for base field).
    pub fn level(&self) -> usize {
        match self {
            TowerElement::Base(_) => 0,
            TowerElement::Extension { level, .. } => *level,
        }
    }

    /// Returns true if this is a base field element.
    pub fn is_base(&self) -> bool {
        matches!(self, TowerElement::Base(_))
    }

    /// Returns the base rational function if this is a base element.
    pub fn as_base(&self) -> Option<&RationalFunction<F>> {
        match self {
            TowerElement::Base(rf) => Some(rf),
            _ => None,
        }
    }

    /// Returns true if this element is zero.
    pub fn is_zero(&self) -> bool {
        match self {
            TowerElement::Base(rf) => rf.is_zero(),
            TowerElement::Extension {
                numerator_coeffs, ..
            } => numerator_coeffs.iter().all(|c| c.is_zero()),
        }
    }

    /// Creates a polynomial element c₀ + c₁θ + c₂θ² + ... at the given level.
    pub fn polynomial(level: usize, coeffs: Vec<TowerElement<F>>) -> Self {
        if coeffs.is_empty() {
            return TowerElement::zero();
        }

        // Check if this is just a constant
        if coeffs.len() == 1 && level == 1 {
            if let Some(base) = coeffs[0].as_base() {
                return TowerElement::Base(base.clone());
            }
        }

        TowerElement::Extension {
            level,
            numerator_coeffs: coeffs,
            denominator_coeffs: vec![TowerElement::one()],
        }
    }

    /// Creates θ (the transcendental) at level 1.
    pub fn theta() -> Self {
        TowerElement::Extension {
            level: 1,
            numerator_coeffs: vec![TowerElement::zero(), TowerElement::one()],
            denominator_coeffs: vec![TowerElement::one()],
        }
    }

    /// Creates θ^n at level 1.
    pub fn theta_pow(n: usize) -> Self {
        if n == 0 {
            return TowerElement::one();
        }

        let mut coeffs = vec![TowerElement::zero(); n + 1];
        coeffs[n] = TowerElement::one();

        TowerElement::polynomial(1, coeffs)
    }
}

/// Operations on tower elements relative to a tower.
impl<F: Field> TowerElement<F> {
    /// Computes the derivative of this element with respect to the tower's derivation.
    ///
    /// For base elements: use standard derivation
    /// For extensions: use the chain rule based on whether θ is log or exp
    pub fn derivative(&self, _tower: &TranscendentalTower) -> TowerElement<F> {
        match self {
            TowerElement::Base(rf) => TowerElement::Base(rf.derivative()),
            TowerElement::Extension { .. } => {
                // TODO: Implement derivative for extension elements
                // This requires knowing the type of θ (log or exp) from the tower
                TowerElement::zero()
            }
        }
    }
}

/// Display for tower elements.
impl<F: Field> std::fmt::Display for TowerElement<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TowerElement::Base(rf) => write!(f, "{}", rf),
            TowerElement::Extension {
                level,
                numerator_coeffs,
                denominator_coeffs,
            } => {
                let theta = format!("θ{}", level);

                // Format numerator
                let mut num_terms = Vec::new();
                for (i, coeff) in numerator_coeffs.iter().enumerate() {
                    if !coeff.is_zero() {
                        let term = match i {
                            0 => format!("({})", coeff),
                            1 => format!("({})*{}", coeff, theta),
                            _ => format!("({})*{}^{}", coeff, theta, i),
                        };
                        num_terms.push(term);
                    }
                }

                let num_str = if num_terms.is_empty() {
                    "0".to_string()
                } else {
                    num_terms.join(" + ")
                };

                // Check if denominator is 1
                let is_denom_one = denominator_coeffs.len() == 1
                    && denominator_coeffs[0]
                        .as_base()
                        .map_or(false, |r| r.is_polynomial() && r.numerator().degree() == 0);

                if is_denom_one {
                    write!(f, "{}", num_str)
                } else {
                    // Format denominator similarly
                    let mut den_terms = Vec::new();
                    for (i, coeff) in denominator_coeffs.iter().enumerate() {
                        if !coeff.is_zero() {
                            let term = match i {
                                0 => format!("({})", coeff),
                                1 => format!("({})*{}", coeff, theta),
                                _ => format!("({})*{}^{}", coeff, theta, i),
                            };
                            den_terms.push(term);
                        }
                    }
                    let den_str = den_terms.join(" + ");

                    write!(f, "({}) / ({})", num_str, den_str)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_base_element() {
        let rf = RationalFunction::new(poly(&[1, 1]), poly(&[1]));
        let elem = TowerElement::from_base(rf);

        assert!(elem.is_base());
        assert_eq!(elem.level(), 0);
    }

    #[test]
    fn test_theta() {
        let theta: TowerElement<Q> = TowerElement::theta();

        assert!(!theta.is_base());
        assert_eq!(theta.level(), 1);
    }

    #[test]
    fn test_theta_power() {
        let theta_sq: TowerElement<Q> = TowerElement::theta_pow(2);

        assert_eq!(theta_sq.level(), 1);
        assert!(!theta_sq.is_zero());
    }

    #[test]
    fn test_polynomial_in_theta() {
        // 1 + θ
        let coeffs = vec![TowerElement::one(), TowerElement::one()];
        let elem: TowerElement<Q> = TowerElement::polynomial(1, coeffs);

        assert_eq!(elem.level(), 1);
    }

    #[test]
    fn test_zero_element() {
        let zero: TowerElement<Q> = TowerElement::zero();
        assert!(zero.is_zero());
    }

    #[test]
    fn test_display() {
        let theta: TowerElement<Q> = TowerElement::theta();
        let s = format!("{}", theta);

        assert!(s.contains("θ1"));
    }
}
