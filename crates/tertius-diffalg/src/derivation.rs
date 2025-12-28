//! Derivation operators for differential algebra.
//!
//! A derivation D on a ring R is a map D: R → R satisfying:
//! - D(a + b) = D(a) + D(b)  (additivity)
//! - D(a · b) = D(a)·b + a·D(b)  (Leibniz rule)
//!
//! For the Risch algorithm, we use the standard derivation d/dx.

use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::traits::Field;

/// A derivation on a differential ring.
pub trait Derivation<T> {
    /// Computes the derivative of an element.
    fn derive(&self, x: &T) -> T;

    /// Checks if an element is a constant (derivative is zero).
    fn is_constant(&self, x: &T) -> bool;
}

/// The standard derivation d/dx on polynomials.
#[derive(Clone, Debug, Default)]
pub struct StandardDerivation;

impl<F: Field> Derivation<DensePoly<F>> for StandardDerivation {
    fn derive(&self, p: &DensePoly<F>) -> DensePoly<F> {
        p.derivative()
    }

    fn is_constant(&self, p: &DensePoly<F>) -> bool {
        p.degree() == 0
    }
}

impl<F: Field> Derivation<RationalFunction<F>> for StandardDerivation {
    fn derive(&self, f: &RationalFunction<F>) -> RationalFunction<F> {
        f.derivative()
    }

    fn is_constant(&self, f: &RationalFunction<F>) -> bool {
        f.is_polynomial() && f.numerator().degree() == 0
    }
}

/// A derivation extended to a transcendental extension.
///
/// For θ = log(u), we have D(θ) = D(u)/u
/// For θ = exp(u), we have D(θ) = D(u)·θ
#[derive(Clone, Debug)]
pub struct ExtendedDerivation<D> {
    /// The base derivation.
    pub base: D,
    /// How to compute D(θ) for each transcendental.
    pub theta_derivatives: Vec<ThetaDerivative>,
}

/// How to compute the derivative of a transcendental θ.
#[derive(Clone, Debug)]
pub enum ThetaDerivative {
    /// θ = log(u): D(θ) = D(u)/u
    Logarithmic {
        /// Index of u in the tower (or representation as polynomial)
        argument_level: usize,
    },
    /// θ = exp(u): D(θ) = D(u)·θ
    Exponential {
        /// Index of u in the tower
        exponent_level: usize,
    },
}

impl<D: Default> Default for ExtendedDerivation<D> {
    fn default() -> Self {
        Self {
            base: D::default(),
            theta_derivatives: Vec::new(),
        }
    }
}

impl<D> ExtendedDerivation<D> {
    /// Creates a new extended derivation.
    pub fn new(base: D) -> Self {
        Self {
            base,
            theta_derivatives: Vec::new(),
        }
    }

    /// Adds a logarithmic extension θ = log(u).
    pub fn add_logarithmic(&mut self, argument_level: usize) {
        self.theta_derivatives
            .push(ThetaDerivative::Logarithmic { argument_level });
    }

    /// Adds an exponential extension θ = exp(u).
    pub fn add_exponential(&mut self, exponent_level: usize) {
        self.theta_derivatives
            .push(ThetaDerivative::Exponential { exponent_level });
    }

    /// Returns the number of transcendental extensions.
    pub fn num_extensions(&self) -> usize {
        self.theta_derivatives.len()
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
    fn test_polynomial_derivation() {
        let d = StandardDerivation;

        // d/dx (x² + 2x + 1) = 2x + 2
        let p = poly(&[1, 2, 1]);
        let dp = d.derive(&p);

        assert_eq!(dp.degree(), 1);
        assert_eq!(dp.coeff(0), q(2));
        assert_eq!(dp.coeff(1), q(2));
    }

    #[test]
    fn test_constant_detection() {
        let d = StandardDerivation;

        let constant = poly(&[5]);
        let non_constant = poly(&[1, 1]);

        assert!(d.is_constant(&constant));
        assert!(!d.is_constant(&non_constant));
    }

    #[test]
    fn test_rational_derivation() {
        let d = StandardDerivation;

        // d/dx (1/x) = -1/x²
        let f = RationalFunction::new(poly(&[1]), poly(&[0, 1]));
        let df = d.derive(&f);

        // Should be -1/x²
        assert_eq!(df.numerator().degree(), 0);
        assert_eq!(df.denominator().degree(), 2);
    }

    #[test]
    fn test_extended_derivation_setup() {
        let mut d = ExtendedDerivation::new(StandardDerivation);

        d.add_logarithmic(0); // θ₁ = log(x)
        d.add_exponential(0); // θ₂ = exp(x)

        assert_eq!(d.num_extensions(), 2);
    }
}
