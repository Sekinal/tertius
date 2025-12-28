//! Error function and related functions.
//!
//! The error function is defined as:
//!
//! erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
//!
//! Related functions:
//! - erfc(x) = 1 - erf(x) (complementary error function)
//! - erfi(x) = -i·erf(ix) (imaginary error function)
//! - Dawson function: D(x) = e^(-x²) ∫₀ˣ e^(t²) dt
//!
//! # Key Properties
//!
//! - erf(0) = 0
//! - erf(∞) = 1
//! - erf(-x) = -erf(x) (odd function)
//! - d/dx erf(x) = (2/√π) e^(-x²)
//!
//! # Integration
//!
//! The canonical non-elementary integral:
//! - ∫ e^(-x²) dx = (√π/2) erf(x)

use tertius_rings::rationals::Q;
use tertius_rings::traits::Field;

/// The error function erf(argument).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ErrorFunction {
    /// Symbolic representation of the argument.
    pub argument: String,
    /// Variant: standard erf, complementary erfc, or imaginary erfi.
    pub variant: ErrorFunctionVariant,
}

/// Variants of the error function.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ErrorFunctionVariant {
    /// Standard error function erf(x)
    Standard,
    /// Complementary error function erfc(x) = 1 - erf(x)
    Complementary,
    /// Imaginary error function erfi(x) = -i·erf(ix)
    Imaginary,
}

impl ErrorFunction {
    /// Creates erf(arg).
    pub fn erf(argument: &str) -> Self {
        Self {
            argument: argument.to_string(),
            variant: ErrorFunctionVariant::Standard,
        }
    }

    /// Creates erfc(arg) = 1 - erf(arg).
    pub fn erfc(argument: &str) -> Self {
        Self {
            argument: argument.to_string(),
            variant: ErrorFunctionVariant::Complementary,
        }
    }

    /// Creates erfi(arg) = -i·erf(i·arg).
    pub fn erfi(argument: &str) -> Self {
        Self {
            argument: argument.to_string(),
            variant: ErrorFunctionVariant::Imaginary,
        }
    }

    /// Returns the derivative of this error function.
    pub fn derivative(&self) -> String {
        match self.variant {
            ErrorFunctionVariant::Standard => {
                format!("(2/√π) * exp(-({}²))", self.argument)
            }
            ErrorFunctionVariant::Complementary => {
                format!("-(2/√π) * exp(-({}²))", self.argument)
            }
            ErrorFunctionVariant::Imaginary => {
                format!("(2/√π) * exp({}²)", self.argument)
            }
        }
    }

    /// Returns the integral ∫erf(x)dx = x·erf(x) + e^(-x²)/√π.
    pub fn integral(&self) -> String {
        match self.variant {
            ErrorFunctionVariant::Standard => {
                format!(
                    "{0} * erf({0}) + exp(-{0}²)/√π",
                    self.argument
                )
            }
            _ => format!("∫ {} dx", self),
        }
    }
}

impl std::fmt::Display for ErrorFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.variant {
            ErrorFunctionVariant::Standard => write!(f, "erf({})", self.argument),
            ErrorFunctionVariant::Complementary => write!(f, "erfc({})", self.argument),
            ErrorFunctionVariant::Imaginary => write!(f, "erfi({})", self.argument),
        }
    }
}

/// Dawson's integral D(x) = e^(-x²) ∫₀ˣ e^(t²) dt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DawsonFunction {
    /// Symbolic representation of the argument.
    pub argument: String,
}

impl DawsonFunction {
    /// Creates D(arg).
    pub fn new(argument: &str) -> Self {
        Self {
            argument: argument.to_string(),
        }
    }

    /// D'(x) = 1 - 2x D(x)
    pub fn derivative(&self) -> String {
        format!("1 - 2*{0}*D({0})", self.argument)
    }
}

impl std::fmt::Display for DawsonFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "D({})", self.argument)
    }
}

/// Computes the error function series expansion.
///
/// erf(x) = (2/√π) Σ_{n=0}^∞ (-1)^n x^(2n+1) / (n! (2n+1))
pub fn error_func_series(x: Q, num_terms: usize) -> Q {
    let mut result = Q::from_integer(0);
    let x_sq = x.clone() * x.clone();
    let mut power = x.clone(); // x^1
    let mut n_factorial = Q::from_integer(1);

    for n in 0..num_terms {
        if n > 0 {
            n_factorial = n_factorial * Q::from_integer(n as i64);
            power = power * x_sq.clone();
        }

        let two_n_plus_1 = Q::from_integer((2 * n + 1) as i64);
        let denom = n_factorial.clone() * two_n_plus_1;

        let sign = if n % 2 == 0 { Q::from_integer(1) } else { Q::from_integer(-1) };
        let term = sign * power.clone() * denom.inv().unwrap();

        result = result + term;
    }

    // Multiply by 2/√π... but we can't represent √π exactly
    // Return the series without the constant factor
    result
}

/// Known integrals involving e^(-x²).
pub struct GaussianIntegrals;

impl GaussianIntegrals {
    /// ∫ e^(-x²) dx = (√π/2) erf(x) + C
    pub fn basic() -> String {
        "(√π/2) * erf(x)".to_string()
    }

    /// ∫ x e^(-x²) dx = -e^(-x²)/2 + C
    pub fn x_times() -> String {
        "-exp(-x²)/2".to_string()
    }

    /// ∫ x² e^(-x²) dx = -(x/2) e^(-x²) + (√π/4) erf(x) + C
    pub fn x_squared_times() -> String {
        "-(x/2)*exp(-x²) + (√π/4)*erf(x)".to_string()
    }

    /// ∫₀^∞ e^(-x²) dx = √π/2
    pub fn definite_zero_to_inf() -> String {
        "√π/2".to_string()
    }

    /// ∫₋∞^∞ e^(-x²) dx = √π (Gaussian integral)
    pub fn gaussian_integral() -> String {
        "√π".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erf_creation() {
        let erf = ErrorFunction::erf("x");
        assert_eq!(format!("{}", erf), "erf(x)");
    }

    #[test]
    fn test_erfc_creation() {
        let erfc = ErrorFunction::erfc("x");
        assert_eq!(format!("{}", erfc), "erfc(x)");
    }

    #[test]
    fn test_erf_derivative() {
        let erf = ErrorFunction::erf("x");
        let deriv = erf.derivative();
        assert!(deriv.contains("exp"));
        assert!(deriv.contains("√π"));
    }

    #[test]
    fn test_dawson_function() {
        let d = DawsonFunction::new("x");
        assert_eq!(format!("{}", d), "D(x)");
    }

    #[test]
    fn test_error_func_series() {
        // erf(0) = 0
        let result = error_func_series(Q::from_integer(0), 10);
        assert_eq!(result, Q::from_integer(0));
    }

    #[test]
    fn test_gaussian_integrals() {
        let basic = GaussianIntegrals::basic();
        assert!(basic.contains("erf"));

        let gaussian = GaussianIntegrals::gaussian_integral();
        assert_eq!(gaussian, "√π");
    }
}
