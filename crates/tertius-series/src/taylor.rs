//! Taylor expansion of functions around a point.
//!
//! Provides the `TaylorExpansion` type for representing Taylor series
//! f(x) = Σᵢ (f⁽ⁱ⁾(a)/i!) * (x-a)ⁱ around point a.

use crate::power_series::{PowerSeries, SeriesCoeff};
use tertius_rings::traits::{Field, Ring};

/// A Taylor expansion around a point.
#[derive(Clone)]
pub struct TaylorExpansion<R: SeriesCoeff> {
    /// The expansion point a.
    pub center: R,
    /// The power series in (x-a).
    pub series: PowerSeries<R>,
}

/// Helper function to create integer n as a ring element.
fn from_usize<R: Ring + Clone>(n: usize) -> R {
    let mut result = R::zero();
    let one = R::one();
    for _ in 0..n {
        result = result + one.clone();
    }
    result
}

impl<R: Field + Clone + Send + Sync + 'static> TaylorExpansion<R> {
    /// Creates a Taylor expansion from derivatives at a point.
    ///
    /// Given f(a), f'(a), f''(a), ..., constructs the Taylor series.
    pub fn from_derivatives(center: R, derivatives: Vec<R>) -> Self {
        let precision = derivatives.len();

        // Compute factorial coefficients
        let mut factorials = vec![R::one()];
        for i in 1..precision {
            let next = factorials[i - 1].clone() * from_usize::<R>(i);
            factorials.push(next);
        }

        // Coefficients are f^(n)(a) / n!
        let coeffs: Vec<R> = derivatives
            .into_iter()
            .enumerate()
            .map(|(i, d)| d.field_div(&factorials[i]))
            .collect();

        Self {
            center,
            series: PowerSeries::from_coeffs(coeffs),
        }
    }

    /// Creates a Taylor expansion from explicit coefficients.
    ///
    /// The coefficients are already divided by n!.
    pub fn from_coefficients(center: R, coefficients: Vec<R>) -> Self {
        Self {
            center,
            series: PowerSeries::from_coeffs(coefficients),
        }
    }

    /// Creates a Taylor expansion from a power series at a point.
    pub fn new(center: R, series: PowerSeries<R>) -> Self {
        Self { center, series }
    }

    /// Returns the expansion center.
    pub fn center(&self) -> &R {
        &self.center
    }

    /// Returns the underlying power series.
    pub fn series(&self) -> &PowerSeries<R> {
        &self.series
    }

    /// Returns the coefficient of (x-a)^n.
    pub fn coeff(&self, n: usize) -> R {
        self.series.coeff(n)
    }

    /// Evaluates the Taylor series at a point x.
    ///
    /// Computes f(x) ≈ Σᵢ aᵢ(x-center)ⁱ
    pub fn evaluate(&self, x: &R) -> R {
        let h = x.clone() - self.center.clone();
        let mut result = R::zero();
        let mut h_power = R::one();

        for i in 0..self.series.precision() {
            result = result + self.series.coeff(i) * h_power.clone();
            h_power = h_power * h.clone();
        }

        result
    }

    /// Returns the precision (number of terms).
    pub fn precision(&self) -> usize {
        self.series.precision()
    }

    /// Shifts the expansion to a new center.
    ///
    /// Given f(x) as Taylor series at a, computes the Taylor series at b.
    pub fn recenter(&self, new_center: R) -> Self {
        // Use binomial theorem to shift
        // If f(x) = Σ aₙ(x-a)ⁿ at center a
        // Then at center b: (x-a) = (x-b) + (b-a)
        // So (x-a)ⁿ = Σₖ C(n,k) (x-b)ᵏ (b-a)^{n-k}

        let precision = self.series.precision();
        let delta = new_center.clone() - self.center.clone();
        let old = self.clone();

        // Compute new coefficients using binomial expansion
        let mut new_coeffs = vec![R::zero(); precision];

        for n in 0..precision {
            let a_n = old.coeff(n);
            if a_n.is_zero() {
                continue;
            }

            for k in 0..=n {
                let binom = binomial::<R>(n, k);
                let delta_pow = pow_ring(&delta, n - k);
                new_coeffs[k] = new_coeffs[k].clone() + binom * a_n.clone() * delta_pow;
            }
        }

        Self {
            center: new_center,
            series: PowerSeries::from_coeffs(new_coeffs),
        }
    }
}

/// Standard Taylor expansions at origin.
impl<R: Field + Clone + Send + Sync + 'static> TaylorExpansion<R> {
    /// Taylor expansion of exp(x) at 0.
    pub fn exp(precision: usize) -> Self {
        Self {
            center: R::zero(),
            series: PowerSeries::exp(precision),
        }
    }

    /// Taylor expansion of sin(x) at 0.
    pub fn sin(precision: usize) -> Self {
        Self {
            center: R::zero(),
            series: PowerSeries::sin(precision),
        }
    }

    /// Taylor expansion of cos(x) at 0.
    pub fn cos(precision: usize) -> Self {
        Self {
            center: R::zero(),
            series: PowerSeries::cos(precision),
        }
    }

    /// Taylor expansion of log(1+x) at 0.
    pub fn log1p(precision: usize) -> Self {
        Self {
            center: R::zero(),
            series: PowerSeries::log1p(precision),
        }
    }

    /// Taylor expansion of 1/(1-x) at 0.
    pub fn geometric(precision: usize) -> Self {
        Self {
            center: R::zero(),
            series: PowerSeries::geometric(precision),
        }
    }
}

/// Computes binomial coefficient C(n, k) as a ring element.
fn binomial<R: Field + Clone>(n: usize, k: usize) -> R {
    if k > n {
        return R::zero();
    }

    // C(n, k) = n! / (k! * (n-k)!)
    // Compute using the formula C(n, k) = C(n, k-1) * (n-k+1) / k

    let mut result = R::one();
    for i in 0..k {
        let numer = from_usize::<R>(n - i);
        let denom = from_usize::<R>(i + 1);
        result = (result * numer).field_div(&denom);
    }
    result
}

/// Computes a^n for a ring element.
fn pow_ring<R: Ring + Clone>(a: &R, n: usize) -> R {
    if n == 0 {
        return R::one();
    }
    let mut result = a.clone();
    for _ in 1..n {
        result = result * a.clone();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64, d: i64) -> Q {
        Q::new(n, d)
    }

    #[test]
    fn test_from_derivatives() {
        // f(x) = exp(x) at x=0: f(0)=1, f'(0)=1, f''(0)=1, ...
        let taylor = TaylorExpansion::from_derivatives(
            q(0, 1),
            vec![q(1, 1), q(1, 1), q(1, 1), q(1, 1)],
        );

        // Coefficients should be 1, 1, 1/2, 1/6
        assert_eq!(taylor.coeff(0), q(1, 1));
        assert_eq!(taylor.coeff(1), q(1, 1));
        assert_eq!(taylor.coeff(2), q(1, 2));
        assert_eq!(taylor.coeff(3), q(1, 6));
    }

    #[test]
    fn test_from_coefficients() {
        let taylor =
            TaylorExpansion::from_coefficients(q(0, 1), vec![q(1, 1), q(2, 1), q(3, 1)]);

        assert_eq!(taylor.coeff(0), q(1, 1));
        assert_eq!(taylor.coeff(1), q(2, 1));
        assert_eq!(taylor.coeff(2), q(3, 1));
    }

    #[test]
    fn test_evaluate() {
        // f(x) = 1 + x + x² at center 0
        let taylor = TaylorExpansion::from_coefficients(
            q(0, 1),
            vec![q(1, 1), q(1, 1), q(1, 1)],
        );

        // f(2) = 1 + 2 + 4 = 7
        assert_eq!(taylor.evaluate(&q(2, 1)), q(7, 1));
    }

    #[test]
    fn test_exp_taylor() {
        let exp: TaylorExpansion<Q> = TaylorExpansion::exp(10);

        // exp(x) = 1 + x + x²/2 + x³/6 + ...
        assert_eq!(exp.coeff(0), q(1, 1));
        assert_eq!(exp.coeff(1), q(1, 1));
        assert_eq!(exp.coeff(2), q(1, 2));
        assert_eq!(exp.coeff(3), q(1, 6));
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial::<Q>(4, 0), q(1, 1));
        assert_eq!(binomial::<Q>(4, 1), q(4, 1));
        assert_eq!(binomial::<Q>(4, 2), q(6, 1));
        assert_eq!(binomial::<Q>(4, 3), q(4, 1));
        assert_eq!(binomial::<Q>(4, 4), q(1, 1));
    }

    #[test]
    fn test_recenter() {
        // f(x) = 1 + x (polynomial) at center 0
        // At center 1: f(x) = f(1 + (x-1)) = 1 + 1 + (x-1) = 2 + (x-1)
        let taylor = TaylorExpansion::from_coefficients(q(0, 1), vec![q(1, 1), q(1, 1)]);

        let recentered = taylor.recenter(q(1, 1));
        assert_eq!(recentered.center(), &q(1, 1));
        assert_eq!(recentered.coeff(0), q(2, 1)); // constant term at new center
        assert_eq!(recentered.coeff(1), q(1, 1)); // linear coefficient unchanged
    }
}
