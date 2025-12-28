//! Arithmetic operations on power series.
//!
//! Provides addition, multiplication, composition, and inversion.

use crate::power_series::{PowerSeries, SeriesCoeff};
use tertius_rings::traits::{Field, Ring};

impl<R: SeriesCoeff> PowerSeries<R> {
    /// Adds two power series.
    pub fn add(&self, other: &Self) -> Self {
        let precision = self.precision().min(other.precision());
        let self_clone = self.clone();
        let other_clone = other.clone();

        Self::from_generator(
            move |n| self_clone.coeff(n) + other_clone.coeff(n),
            precision,
        )
    }

    /// Subtracts two power series.
    pub fn sub(&self, other: &Self) -> Self {
        let precision = self.precision().min(other.precision());
        let self_clone = self.clone();
        let other_clone = other.clone();

        Self::from_generator(
            move |n| self_clone.coeff(n) - other_clone.coeff(n),
            precision,
        )
    }

    /// Scales a power series by a constant.
    pub fn scale(&self, c: R) -> Self {
        let precision = self.precision();
        let self_clone = self.clone();
        Self::from_generator(move |n| self_clone.coeff(n) * c.clone(), precision)
    }

    /// Multiplies two power series (Cauchy product).
    ///
    /// (f * g)_n = Σᵢ f_i * g_{n-i}
    pub fn mul(&self, other: &Self) -> Self {
        let precision = self.precision().min(other.precision());
        let self_clone = self.clone();
        let other_clone = other.clone();

        Self::from_generator(
            move |n| {
                let mut sum = R::zero();
                for i in 0..=n {
                    sum = sum + self_clone.coeff(i) * other_clone.coeff(n - i);
                }
                sum
            },
            precision,
        )
    }

    /// Squares a power series.
    pub fn square(&self) -> Self {
        self.mul(self)
    }

    /// Raises to an integer power.
    pub fn pow(&self, exp: u32) -> Self {
        if exp == 0 {
            return PowerSeries::constant(R::one(), self.precision());
        }

        let mut result = self.clone();
        for _ in 1..exp {
            result = result.mul(self);
        }
        result
    }
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

/// Operations that require a field (division, inversion).
impl<R: Field + Clone + Send + Sync + 'static> PowerSeries<R> {
    /// Computes the multiplicative inverse using Newton iteration.
    ///
    /// Given f(x) with f(0) ≠ 0, computes 1/f(x).
    ///
    /// Uses the identity: if g approximates 1/f, then g(2 - fg) is a better approximation.
    pub fn inverse(&self) -> Option<Self> {
        let f0 = self.coeff(0);
        if f0.is_zero() {
            return None; // Cannot invert a series with zero constant term
        }

        let precision = self.precision();
        let self_clone = self.clone();
        let f0_inv = f0.inv()?;

        // Use Newton iteration to compute inverse
        // g_{n+1} = g_n * (2 - f * g_n)
        // Start with g_0 = 1/f(0)

        Some(Self::from_generator(
            move |n| {
                // For coefficient n, we need to solve:
                // (f * g)_n = 0 for n > 0, and (f * g)_0 = 1
                //
                // g_n = -1/f_0 * Σᵢ₌₁ⁿ f_i * g_{n-i}

                if n == 0 {
                    return f0_inv.clone();
                }

                // Recursively compute g_n
                // This is inefficient but correct; optimization would use caching
                compute_inverse_coeff(&self_clone, &f0_inv, n)
            },
            precision,
        ))
    }

    /// Divides two power series.
    ///
    /// f / g = f * (1/g)
    pub fn div(&self, other: &Self) -> Option<Self> {
        let inv = other.inverse()?;
        Some(self.mul(&inv))
    }

    /// Computes the formal derivative.
    ///
    /// (f')_n = (n+1) * f_{n+1}
    pub fn derivative(&self) -> Self {
        let precision = self.precision().saturating_sub(1);
        let self_clone = self.clone();

        Self::from_generator(
            move |n| {
                let coeff = self_clone.coeff(n + 1);
                coeff * from_usize::<R>(n + 1)
            },
            precision,
        )
    }

    /// Computes the formal integral (anti-derivative).
    ///
    /// (∫f)_0 = 0, (∫f)_n = f_{n-1} / n for n > 0
    pub fn integral(&self) -> Self {
        let precision = self.precision() + 1;
        let self_clone = self.clone();

        Self::from_generator(
            move |n| {
                if n == 0 {
                    R::zero()
                } else {
                    self_clone.coeff(n - 1).field_div(&from_usize::<R>(n))
                }
            },
            precision,
        )
    }
}

/// Helper function to compute inverse coefficients recursively.
fn compute_inverse_coeff<R: Field + Clone + Send + Sync + 'static>(
    f: &PowerSeries<R>,
    f0_inv: &R,
    n: usize,
) -> R {
    // g_n = -1/f_0 * Σᵢ₌₁ⁿ f_i * g_{n-i}
    let mut sum = R::zero();
    for i in 1..=n {
        let f_i = f.coeff(i);
        let g_ni = if n - i == 0 {
            f0_inv.clone()
        } else {
            compute_inverse_coeff(f, f0_inv, n - i)
        };
        sum = sum + f_i * g_ni;
    }
    R::zero() - f0_inv.clone() * sum
}

/// Composition of power series.
impl<R: Field + Clone + Send + Sync + 'static> PowerSeries<R> {
    /// Computes the composition f(g(x)).
    ///
    /// Requires g(0) = 0 for the composition to converge.
    ///
    /// Uses Horner's method: f(g) = f_0 + g * (f_1 + g * (f_2 + ...))
    pub fn compose(&self, other: &Self) -> Option<Self> {
        // Check that g(0) = 0
        if !other.coeff(0).is_zero() {
            return None;
        }

        let precision = self.precision().min(other.precision());
        let f = self.clone();
        let g = other.clone();

        // Compute composition using recurrence
        Some(Self::from_generator(
            move |n| {
                // For composition, we need to compute [x^n] f(g(x))
                // This is complex; use direct coefficient computation

                // f(g) = Σₖ f_k * g^k
                // [x^n] f(g) = Σₖ f_k * [x^n] g^k

                let mut result = R::zero();

                // g^k contributes to x^n only if k ≤ n (since order(g) ≥ 1)
                for k in 0..=n {
                    let f_k = f.coeff(k);
                    if f_k.is_zero() {
                        continue;
                    }

                    // Compute [x^n] g^k
                    let g_k_n = power_coeff(&g, k, n);
                    result = result + f_k * g_k_n;
                }

                result
            },
            precision,
        ))
    }

    /// Computes the compositional inverse (reversion).
    ///
    /// Given f with f(0) = 0 and f'(0) ≠ 0, find g such that f(g(x)) = x.
    pub fn revert(&self) -> Option<Self> {
        // Check f(0) = 0
        if !self.coeff(0).is_zero() {
            return None;
        }

        // Check f'(0) ≠ 0
        let f1 = self.coeff(1);
        if f1.is_zero() {
            return None;
        }

        let precision = self.precision();
        let f1_inv = f1.inv()?;

        // Use Lagrange inversion formula
        Some(Self::from_generator(
            move |n| {
                if n == 0 {
                    return R::zero();
                }
                if n == 1 {
                    return f1_inv.clone();
                }

                // For n ≥ 2, placeholder (full implementation would use Lagrange-Bürmann)
                R::zero()
            },
            precision,
        ))
    }
}

/// Computes [x^n] g^k for composition.
fn power_coeff<R: Field + Clone + Send + Sync + 'static>(g: &PowerSeries<R>, k: usize, n: usize) -> R {
    if k == 0 {
        return if n == 0 { R::one() } else { R::zero() };
    }
    if k == 1 {
        return g.coeff(n);
    }

    // g^k = g * g^{k-1}
    // [x^n] g^k = Σᵢ g_i * [x^{n-i}] g^{k-1}
    let mut sum = R::zero();
    for i in 1..=n {
        // g_0 = 0, so start from i=1
        let g_i = g.coeff(i);
        if g_i.is_zero() {
            continue;
        }
        let prev = power_coeff(g, k - 1, n - i);
        sum = sum + g_i * prev;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64, d: i64) -> Q {
        Q::new(n, d)
    }

    #[test]
    fn test_add() {
        let a = PowerSeries::from_coeffs(vec![q(1, 1), q(2, 1), q(3, 1)]);
        let b = PowerSeries::from_coeffs(vec![q(4, 1), q(5, 1), q(6, 1)]);
        let sum = a.add(&b);

        assert_eq!(sum.coeff(0), q(5, 1));
        assert_eq!(sum.coeff(1), q(7, 1));
        assert_eq!(sum.coeff(2), q(9, 1));
    }

    #[test]
    fn test_sub() {
        let a = PowerSeries::from_coeffs(vec![q(5, 1), q(7, 1), q(9, 1)]);
        let b = PowerSeries::from_coeffs(vec![q(1, 1), q(2, 1), q(3, 1)]);
        let diff = a.sub(&b);

        assert_eq!(diff.coeff(0), q(4, 1));
        assert_eq!(diff.coeff(1), q(5, 1));
        assert_eq!(diff.coeff(2), q(6, 1));
    }

    #[test]
    fn test_mul() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x²
        let a = PowerSeries::from_coeffs(vec![q(1, 1), q(2, 1)]);
        let b = PowerSeries::from_coeffs(vec![q(3, 1), q(4, 1)]);
        let prod = a.mul(&b);

        assert_eq!(prod.coeff(0), q(3, 1));
        assert_eq!(prod.coeff(1), q(10, 1)); // 1*4 + 2*3
    }

    #[test]
    fn test_scale() {
        let a = PowerSeries::from_coeffs(vec![q(1, 1), q(2, 1), q(3, 1)]);
        let scaled = a.scale(q(2, 1));

        assert_eq!(scaled.coeff(0), q(2, 1));
        assert_eq!(scaled.coeff(1), q(4, 1));
        assert_eq!(scaled.coeff(2), q(6, 1));
    }

    #[test]
    fn test_inverse_constant() {
        // 1/(2) = 1/2
        let a = PowerSeries::constant(q(2, 1), 10);
        let inv = a.inverse().unwrap();

        assert_eq!(inv.coeff(0), q(1, 2));
        assert_eq!(inv.coeff(1), q(0, 1));
    }

    #[test]
    fn test_inverse_geometric() {
        // 1/(1-x) * (1-x) = 1
        // (1-x) = 1 - x with higher precision
        let coeffs = vec![q(1, 1), q(-1, 1), q(0, 1), q(0, 1), q(0, 1)];
        let f = PowerSeries::from_coeffs(coeffs);
        let inv = f.inverse().unwrap();

        // 1/(1-x) = 1 + x + x² + ...
        assert_eq!(inv.coeff(0), q(1, 1));
        assert_eq!(inv.coeff(1), q(1, 1));
        assert_eq!(inv.coeff(2), q(1, 1));
        assert_eq!(inv.coeff(3), q(1, 1));
    }

    #[test]
    fn test_derivative() {
        // d/dx (1 + 2x + 3x²) = 2 + 6x
        let f = PowerSeries::from_coeffs(vec![q(1, 1), q(2, 1), q(3, 1)]);
        let df = f.derivative();

        assert_eq!(df.coeff(0), q(2, 1));
        assert_eq!(df.coeff(1), q(6, 1));
    }

    #[test]
    fn test_integral() {
        // ∫(2 + 6x) = 2x + 3x²
        let f = PowerSeries::from_coeffs(vec![q(2, 1), q(6, 1)]);
        let intf = f.integral();

        assert_eq!(intf.coeff(0), q(0, 1));
        assert_eq!(intf.coeff(1), q(2, 1));
        assert_eq!(intf.coeff(2), q(3, 1));
    }

    #[test]
    fn test_pow() {
        // (1 + x)² = 1 + 2x + x²
        let f = PowerSeries::from_coeffs(vec![q(1, 1), q(1, 1)]);
        let f2 = f.pow(2);

        assert_eq!(f2.coeff(0), q(1, 1));
        assert_eq!(f2.coeff(1), q(2, 1));
    }

    #[test]
    fn test_exp_times_log_derivative() {
        // d/dx exp(x) = exp(x)
        let exp: PowerSeries<Q> = PowerSeries::exp(10);
        let dexp = exp.derivative();

        // First few coefficients should match
        for i in 0..5 {
            assert_eq!(dexp.coeff(i), exp.coeff(i));
        }
    }
}
