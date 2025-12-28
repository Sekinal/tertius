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

/// Arithmetic operations on tower elements.
impl<F: Field> std::ops::Add for TowerElement<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (TowerElement::Base(a), TowerElement::Base(b)) => TowerElement::Base(a + b),

            // Zero + Extension = Extension
            (a, b) if a.is_zero() => b,
            (a, b) if b.is_zero() => a,

            // Base + Extension: treat base as constant term
            (TowerElement::Base(scalar), TowerElement::Extension { level, mut numerator_coeffs, denominator_coeffs }) |
            (TowerElement::Extension { level, mut numerator_coeffs, denominator_coeffs }, TowerElement::Base(scalar)) => {
                // Add scalar to the constant term (coefficient of θ⁰)
                let is_denom_one = denominator_coeffs.len() == 1 && denominator_coeffs[0].is_one_like();
                if is_denom_one {
                    if numerator_coeffs.is_empty() {
                        numerator_coeffs.push(TowerElement::Base(scalar));
                    } else {
                        numerator_coeffs[0] = numerator_coeffs[0].clone() + TowerElement::Base(scalar);
                    }

                    // Trim trailing zeros
                    while numerator_coeffs.last().map_or(false, |c| c.is_zero()) {
                        numerator_coeffs.pop();
                    }

                    if numerator_coeffs.is_empty() {
                        TowerElement::zero()
                    } else {
                        TowerElement::Extension {
                            level,
                            numerator_coeffs,
                            denominator_coeffs,
                        }
                    }
                } else {
                    TowerElement::zero()
                }
            },

            (TowerElement::Extension { level: l1, numerator_coeffs: n1, denominator_coeffs: d1 },
             TowerElement::Extension { level: l2, numerator_coeffs: n2, denominator_coeffs: d2 }) if l1 == l2 => {
                // (n1/d1) + (n2/d2) = (n1*d2 + n2*d1) / (d1*d2)
                // For simplicity, assume both have denominator 1 for now
                let is_d1_one = d1.len() == 1 && d1[0].is_one_like();
                let is_d2_one = d2.len() == 1 && d2[0].is_one_like();

                if is_d1_one && is_d2_one {
                    // Simple polynomial addition
                    let max_len = n1.len().max(n2.len());
                    let mut result = Vec::with_capacity(max_len);

                    for i in 0..max_len {
                        let c1 = n1.get(i).cloned().unwrap_or_else(TowerElement::zero);
                        let c2 = n2.get(i).cloned().unwrap_or_else(TowerElement::zero);
                        result.push(c1 + c2);
                    }

                    // Trim trailing zeros
                    while result.last().map_or(false, |c| c.is_zero()) {
                        result.pop();
                    }

                    if result.is_empty() {
                        TowerElement::zero()
                    } else {
                        TowerElement::Extension {
                            level: l1,
                            numerator_coeffs: result,
                            denominator_coeffs: vec![TowerElement::one()],
                        }
                    }
                } else {
                    // Full rational function addition - not implemented yet
                    TowerElement::zero()
                }
            },
            // Different levels or other cases - not implemented
            _ => TowerElement::zero(),
        }
    }
}

impl<F: Field> std::ops::Mul for TowerElement<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (TowerElement::Base(a), TowerElement::Base(b)) => TowerElement::Base(a * b),
            (TowerElement::Extension { level: l1, numerator_coeffs: n1, denominator_coeffs: d1 },
             TowerElement::Extension { level: l2, numerator_coeffs: n2, denominator_coeffs: d2 }) if l1 == l2 => {
                // Polynomial multiplication (assuming denominator = 1)
                let is_d1_one = d1.len() == 1 && d1[0].is_one_like();
                let is_d2_one = d2.len() == 1 && d2[0].is_one_like();

                if is_d1_one && is_d2_one && !n1.is_empty() && !n2.is_empty() {
                    let result_len = n1.len() + n2.len() - 1;
                    let mut result = vec![TowerElement::zero(); result_len];

                    for (i, c1) in n1.iter().enumerate() {
                        for (j, c2) in n2.iter().enumerate() {
                            let prod = c1.clone() * c2.clone();
                            result[i + j] = result[i + j].clone() + prod;
                        }
                    }

                    // Trim trailing zeros
                    while result.last().map_or(false, |c| c.is_zero()) {
                        result.pop();
                    }

                    if result.is_empty() {
                        TowerElement::zero()
                    } else {
                        TowerElement::Extension {
                            level: l1,
                            numerator_coeffs: result,
                            denominator_coeffs: vec![TowerElement::one()],
                        }
                    }
                } else {
                    TowerElement::zero()
                }
            },
            // Base * Extension: multiply each coefficient
            (TowerElement::Base(scalar), TowerElement::Extension { level, numerator_coeffs, denominator_coeffs }) |
            (TowerElement::Extension { level, numerator_coeffs, denominator_coeffs }, TowerElement::Base(scalar)) => {
                let new_num: Vec<_> = numerator_coeffs.into_iter()
                    .map(|c| c * TowerElement::Base(scalar.clone()))
                    .collect();
                TowerElement::Extension {
                    level,
                    numerator_coeffs: new_num,
                    denominator_coeffs,
                }
            },
            _ => TowerElement::zero(),
        }
    }
}

impl<F: Field> TowerElement<F> {
    /// Checks if this element is "like one" (for denominator checks).
    fn is_one_like(&self) -> bool {
        match self {
            TowerElement::Base(rf) => {
                // Check if rf = 1/1 (polynomial with value 1)
                rf.is_polynomial()
                    && rf.numerator().degree() == 0
                    && rf.numerator().leading_coeff().is_one()
            },
            _ => false,
        }
    }

    /// Scalar multiplication by an integer.
    pub fn scale(&self, n: i64) -> Self {
        match self {
            TowerElement::Base(rf) => {
                // rf * n
                let scalar = RationalFunction::from_poly(DensePoly::new(vec![F::one().mul_by_scalar(n)]));
                TowerElement::Base(rf.clone() * scalar)
            },
            TowerElement::Extension { level, numerator_coeffs, denominator_coeffs } => {
                TowerElement::Extension {
                    level: *level,
                    numerator_coeffs: numerator_coeffs.iter().map(|c| c.scale(n)).collect(),
                    denominator_coeffs: denominator_coeffs.clone(),
                }
            }
        }
    }
}

/// Operations on tower elements relative to a tower.
impl<F: Field> TowerElement<F> {
    /// Computes the derivative of this element with respect to the tower's derivation.
    ///
    /// For base elements: use standard derivation d/dx
    /// For extensions: use the chain rule based on whether θ is log or exp
    ///
    /// For θ = log(u): θ' = u'/u
    /// For θ = exp(u): θ' = u'·θ
    ///
    /// For polynomial p(θ) = Σ cᵢ θⁱ:
    /// p'(θ) = Σ (cᵢ' θⁱ + i·cᵢ·θⁱ⁻¹·θ')
    pub fn derivative(&self, tower: &TranscendentalTower) -> TowerElement<F> {
        match self {
            TowerElement::Base(rf) => TowerElement::Base(rf.derivative()),
            TowerElement::Extension { level, numerator_coeffs, denominator_coeffs } => {
                // Get the extension type from the tower
                let tower_level = match tower.level(*level - 1) {
                    Some(l) => l,
                    None => return TowerElement::zero(),
                };

                // Check if denominator is 1 (polynomial case)
                let is_polynomial = denominator_coeffs.len() == 1
                    && denominator_coeffs[0].is_one_like();

                if !is_polynomial {
                    // Rational function case: use quotient rule
                    // (p/q)' = (p'q - pq')/q²
                    // This requires more complex arithmetic - return zero for now
                    return TowerElement::zero();
                }

                // Polynomial case: p(θ) = Σ cᵢ θⁱ
                // p'(θ) = Σ (cᵢ' θⁱ + i·cᵢ·θⁱ⁻¹·θ')
                //       = Σ cᵢ' θⁱ + θ' · Σ i·cᵢ·θⁱ⁻¹

                // Compute θ' based on extension type
                let theta_derivative = self.compute_theta_derivative(tower_level, *level, tower);

                // Part 1: Σ cᵢ' θⁱ (coefficient derivatives)
                let coeff_deriv_part = self.derivative_coeff_part(numerator_coeffs, *level, tower);

                // Part 2: θ' · Σ i·cᵢ·θⁱ⁻¹ (chain rule part)
                let chain_rule_part = self.derivative_chain_part(numerator_coeffs, *level, &theta_derivative, tower);

                // Sum the two parts
                coeff_deriv_part + chain_rule_part
            }
        }
    }

    /// Computes θ' for the given extension level.
    fn compute_theta_derivative(&self, tower_level: &crate::tower::TowerLevel, level: usize, _tower: &TranscendentalTower) -> TowerElement<F> {
        match &tower_level.extension_type {
            crate::tower::TranscendentalType::Logarithmic { argument_repr, .. } => {
                // θ = log(u), θ' = u'/u
                // For θ = log(x): θ' = 1/x
                if argument_repr == "x" {
                    // 1/x as a base rational function
                    let one = DensePoly::new(vec![F::one()]);
                    let x = DensePoly::new(vec![F::zero(), F::one()]);
                    TowerElement::Base(RationalFunction::new(one, x))
                } else {
                    // For other arguments, we'd need to parse and compute
                    // For now, return a placeholder (zero means derivative fails gracefully)
                    TowerElement::zero()
                }
            },
            crate::tower::TranscendentalType::Exponential { exponent_repr, .. } => {
                // θ = exp(u), θ' = u'·θ
                // For θ = exp(x): θ' = 1·θ = θ
                if exponent_repr == "x" {
                    // θ itself
                    TowerElement::Extension {
                        level,
                        numerator_coeffs: vec![TowerElement::zero(), TowerElement::one()],
                        denominator_coeffs: vec![TowerElement::one()],
                    }
                } else {
                    // Would need to compute u' for general case
                    TowerElement::zero()
                }
            }
        }
    }

    /// Computes Σ cᵢ' θⁱ (the coefficient derivative part).
    fn derivative_coeff_part(&self, coeffs: &[TowerElement<F>], level: usize, tower: &TranscendentalTower) -> TowerElement<F> {
        let mut result_coeffs: Vec<TowerElement<F>> = Vec::new();

        for c in coeffs {
            result_coeffs.push(c.derivative(tower));
        }

        // Trim trailing zeros
        while result_coeffs.last().map_or(false, |c| c.is_zero()) {
            result_coeffs.pop();
        }

        if result_coeffs.is_empty() {
            TowerElement::zero()
        } else {
            TowerElement::Extension {
                level,
                numerator_coeffs: result_coeffs,
                denominator_coeffs: vec![TowerElement::one()],
            }
        }
    }

    /// Computes θ' · Σ i·cᵢ·θⁱ⁻¹ (the chain rule part).
    fn derivative_chain_part(&self, coeffs: &[TowerElement<F>], level: usize, theta_deriv: &TowerElement<F>, _tower: &TranscendentalTower) -> TowerElement<F> {
        if theta_deriv.is_zero() {
            return TowerElement::zero();
        }

        // Compute p'(θ) = Σ i·cᵢ·θⁱ⁻¹ (formal derivative of polynomial)
        let mut formal_deriv_coeffs: Vec<TowerElement<F>> = Vec::new();

        for (i, c) in coeffs.iter().enumerate().skip(1) {
            // i * c_{i} becomes coefficient of θ^{i-1}
            formal_deriv_coeffs.push(c.scale(i as i64));
        }

        // Trim trailing zeros
        while formal_deriv_coeffs.last().map_or(false, |c| c.is_zero()) {
            formal_deriv_coeffs.pop();
        }

        if formal_deriv_coeffs.is_empty() {
            return TowerElement::zero();
        }

        let formal_deriv = TowerElement::Extension {
            level,
            numerator_coeffs: formal_deriv_coeffs,
            denominator_coeffs: vec![TowerElement::one()],
        };

        // Return θ' * formal_derivative
        theta_deriv.clone() * formal_deriv
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

    #[test]
    fn test_base_derivative() {
        // d/dx(x^2) = 2x
        let tower = TranscendentalTower::new("x");
        let x_squared = TowerElement::from_base(RationalFunction::from_poly(poly(&[0, 0, 1])));

        let deriv = x_squared.derivative(&tower);

        if let TowerElement::Base(rf) = deriv {
            // Should be 2x
            assert_eq!(rf.numerator().degree(), 1);
            assert_eq!(rf.numerator().coeff(1), q(2));
        } else {
            panic!("Expected base element");
        }
    }

    #[test]
    fn test_log_extension_derivative() {
        // θ = log(x), so θ' = 1/x
        // d/dx(θ²) = 2θ · θ' = 2θ/x
        let tower = TranscendentalTower::single_log("x");
        let theta_squared: TowerElement<Q> = TowerElement::theta_pow(2);

        let deriv = theta_squared.derivative(&tower);

        // Should be 2θ · (1/x) = (2/x) · θ
        // This is an extension with coefficient 2/x at θ¹
        assert!(!deriv.is_zero());
        assert!(!deriv.is_base());
    }

    #[test]
    fn test_exp_extension_derivative() {
        // θ = exp(x), so θ' = θ
        // d/dx(θ²) = 2θ · θ' = 2θ · θ = 2θ²
        let tower = TranscendentalTower::single_exp("x");
        let theta_squared: TowerElement<Q> = TowerElement::theta_pow(2);

        let deriv = theta_squared.derivative(&tower);

        // Should be 2θ²
        if let TowerElement::Extension { numerator_coeffs, .. } = deriv {
            // Check coefficient at θ² is 2
            assert!(numerator_coeffs.len() >= 3);
            if let TowerElement::Base(rf) = &numerator_coeffs[2] {
                assert_eq!(rf.numerator().coeff(0), q(2));
            }
        } else {
            panic!("Expected extension element");
        }
    }

    #[test]
    fn test_addition() {
        // θ + θ = 2θ
        let theta: TowerElement<Q> = TowerElement::theta();
        let sum = theta.clone() + theta;

        if let TowerElement::Extension { numerator_coeffs, .. } = sum {
            // Coefficient at θ¹ should be 2
            assert!(numerator_coeffs.len() >= 2);
            if let TowerElement::Base(rf) = &numerator_coeffs[1] {
                assert_eq!(rf.numerator().coeff(0), q(2));
            }
        } else {
            panic!("Expected extension");
        }
    }

    #[test]
    fn test_multiplication() {
        // θ * θ = θ²
        let theta: TowerElement<Q> = TowerElement::theta();
        let prod = theta.clone() * theta;

        if let TowerElement::Extension { numerator_coeffs, .. } = prod {
            // Should have coefficient 1 at θ²
            assert_eq!(numerator_coeffs.len(), 3); // [0, 0, 1]
            assert!(numerator_coeffs[0].is_zero());
            assert!(numerator_coeffs[1].is_zero());
            assert!(numerator_coeffs[2].is_one_like());
        } else {
            panic!("Expected extension");
        }
    }

    #[test]
    fn test_scale() {
        // 3 * θ
        let theta: TowerElement<Q> = TowerElement::theta();
        let scaled = theta.scale(3);

        if let TowerElement::Extension { numerator_coeffs, .. } = scaled {
            // Coefficient at θ¹ should be 3
            if let TowerElement::Base(rf) = &numerator_coeffs[1] {
                assert_eq!(rf.numerator().coeff(0), q(3));
            }
        } else {
            panic!("Expected extension");
        }
    }
}
