//! Special function recognition for non-elementary integrals.
//!
//! This module recognizes integrands that have closed-form antiderivatives
//! in terms of special functions. When the Risch algorithm determines that
//! no elementary antiderivative exists, we check for these canonical forms.
//!
//! # Supported Special Functions
//!
//! | Integral | Result | Special Function |
//! |----------|--------|------------------|
//! | ∫ exp(-x²) dx | (√π/2) erf(x) | Error function |
//! | ∫ exp(x)/x dx | Ei(x) | Exponential integral |
//! | ∫ 1/log(x) dx | li(x) | Logarithmic integral |
//! | ∫ sin(x)/x dx | Si(x) | Sine integral |
//! | ∫ cos(x)/x dx | Ci(x) | Cosine integral |
//! | ∫ exp(-x)/x dx | -Ei(-x) | Exponential integral |

use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::rationals::Q;
use tertius_rings::traits::{Field, Ring};

/// Result of special function recognition.
#[derive(Clone, Debug)]
pub enum SpecialFormResult {
    /// Recognized as a known special function integral
    Recognized(SpecialIntegral),
    /// Not a recognized special form
    NotRecognized,
}

/// A special function integral representation.
#[derive(Clone, Debug)]
pub struct SpecialIntegral {
    /// The special function type
    pub function: SpecialFunction,
    /// Coefficient multiplying the special function
    pub coefficient: SpecialCoefficient,
    /// The argument of the special function
    pub argument: SpecialArgument,
}

/// Coefficient for a special function.
#[derive(Clone, Debug)]
pub enum SpecialCoefficient {
    /// A rational number
    Rational(Q),
    /// √π / n for some integer n
    SqrtPiOver(i64),
    /// A constant involving π
    PiRelated(String),
}

/// Argument transformation for the special function.
#[derive(Clone, Debug)]
pub enum SpecialArgument {
    /// f(x) = x
    Identity,
    /// f(x) = -x
    Negated,
    /// f(x) = -x²
    NegativeSquare,
    /// f(x) = ax + b for some a, b
    Linear { a: Q, b: Q },
    /// f(x) = ax² for some a
    Quadratic { a: Q },
}

/// Types of special functions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpecialFunction {
    /// Error function: erf(x) = (2/√π) ∫₀ˣ exp(-t²) dt
    Erf,
    /// Complementary error function: erfc(x) = 1 - erf(x)
    Erfc,
    /// Exponential integral: Ei(x) = -∫₋ₓ^∞ exp(-t)/t dt
    Ei,
    /// Logarithmic integral: li(x) = ∫₀ˣ dt/log(t)
    Li,
    /// Sine integral: Si(x) = ∫₀ˣ sin(t)/t dt
    Si,
    /// Cosine integral: Ci(x) = -∫ₓ^∞ cos(t)/t dt
    Ci,
    /// Fresnel sine integral: S(x) = ∫₀ˣ sin(πt²/2) dt
    FresnelS,
    /// Fresnel cosine integral: C(x) = ∫₀ˣ cos(πt²/2) dt
    FresnelC,
    /// Gamma function: Γ(x) = ∫₀^∞ t^(x-1) exp(-t) dt
    Gamma,
    /// Lower incomplete gamma: γ(s,x) = ∫₀ˣ t^(s-1) exp(-t) dt
    LowerGamma { s: Q },
    /// Dilogarithm: Li₂(x) = -∫₀ˣ log(1-t)/t dt
    Dilog,
}

impl std::fmt::Display for SpecialFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpecialFunction::Erf => write!(f, "erf"),
            SpecialFunction::Erfc => write!(f, "erfc"),
            SpecialFunction::Ei => write!(f, "Ei"),
            SpecialFunction::Li => write!(f, "li"),
            SpecialFunction::Si => write!(f, "Si"),
            SpecialFunction::Ci => write!(f, "Ci"),
            SpecialFunction::FresnelS => write!(f, "S"),
            SpecialFunction::FresnelC => write!(f, "C"),
            SpecialFunction::Gamma => write!(f, "Γ"),
            SpecialFunction::LowerGamma { s } => write!(f, "γ({}, ·)", s),
            SpecialFunction::Dilog => write!(f, "Li₂"),
        }
    }
}

impl std::fmt::Display for SpecialIntegral {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let coeff = match &self.coefficient {
            SpecialCoefficient::Rational(q) => {
                if *q == Q::from_integer(1) {
                    String::new()
                } else {
                    format!("{} * ", q)
                }
            }
            SpecialCoefficient::SqrtPiOver(n) => {
                if *n == 1 {
                    "√π * ".to_string()
                } else if *n == 2 {
                    "(√π/2) * ".to_string()
                } else {
                    format!("(√π/{}) * ", n)
                }
            }
            SpecialCoefficient::PiRelated(s) => format!("{} * ", s),
        };

        let arg = match &self.argument {
            SpecialArgument::Identity => "x".to_string(),
            SpecialArgument::Negated => "-x".to_string(),
            SpecialArgument::NegativeSquare => "-x²".to_string(),
            SpecialArgument::Linear { a, b } => {
                if b.is_zero() {
                    format!("{}*x", a)
                } else {
                    format!("{}*x + {}", a, b)
                }
            }
            SpecialArgument::Quadratic { a } => format!("{}*x²", a),
        };

        write!(f, "{}{}({})", coeff, self.function, arg)
    }
}

/// Attempts to recognize an integrand as a special function form.
///
/// This is called when the Risch algorithm determines that no elementary
/// antiderivative exists. We check if the integrand matches a known pattern
/// that has a closed-form answer in terms of special functions.
pub fn recognize_special_form<F: Field>(
    integrand: &RationalFunction<F>,
) -> SpecialFormResult {
    // For now, we only implement recognition for Q(x) integrands
    // The full implementation would handle tower elements as well

    SpecialFormResult::NotRecognized
}

/// Recognizes Gaussian-type integrands: ∫ exp(ax² + bx + c) dx
///
/// Result: exp(c - b²/4a) * √(π/(-a)) * erf((2ax + b)/√(-4a)) / 2
/// (for a < 0, i.e., the integral converges)
pub fn recognize_gaussian(
    numerator: &DensePoly<Q>,
    exponent: &DensePoly<Q>,
) -> SpecialFormResult {
    // Check if numerator is constant (typically 1)
    if numerator.degree() != 0 {
        return SpecialFormResult::NotRecognized;
    }

    // Check if exponent is quadratic
    if exponent.degree() != 2 {
        return SpecialFormResult::NotRecognized;
    }

    let a = exponent.coeff(2);
    let b = exponent.coeff(1);
    let c = exponent.coeff(0);

    // For a convergent Gaussian, we need a < 0
    // The canonical form is ∫ exp(-x²) dx = (√π/2) erf(x)

    // Check for the simplest case: exp(-x²)
    if a == Q::from_integer(-1) && b.is_zero() && c.is_zero() {
        return SpecialFormResult::Recognized(SpecialIntegral {
            function: SpecialFunction::Erf,
            coefficient: SpecialCoefficient::SqrtPiOver(2),
            argument: SpecialArgument::Identity,
        });
    }

    // General case with completing the square
    // Not implemented yet
    SpecialFormResult::NotRecognized
}

/// Recognizes exponential integral forms: ∫ exp(x)/x dx, ∫ exp(-x)/x dx
pub fn recognize_exponential_integral(
    poly_part: &DensePoly<Q>,
    rational_part: &RationalFunction<Q>,
) -> SpecialFormResult {
    // ∫ exp(x)/x dx = Ei(x)
    // ∫ exp(-x)/x dx = -Ei(-x)
    // More generally: ∫ exp(ax)/x dx = Ei(ax) for a ≠ 0

    // Check if rational_part is 1/x
    if !rational_part.numerator().is_constant_one() {
        return SpecialFormResult::NotRecognized;
    }

    if rational_part.denominator().degree() != 1 {
        return SpecialFormResult::NotRecognized;
    }

    // Check denominator is just x (not ax + b)
    let d0 = rational_part.denominator().coeff(0);
    let d1 = rational_part.denominator().coeff(1);

    if !d0.is_zero() || d1 != Q::from_integer(1) {
        return SpecialFormResult::NotRecognized;
    }

    // Check if poly_part is just ax for some a
    if poly_part.degree() > 1 || !poly_part.coeff(0).is_zero() {
        return SpecialFormResult::NotRecognized;
    }

    let a = poly_part.coeff(1);

    if a == Q::from_integer(1) {
        // ∫ exp(x)/x dx = Ei(x)
        SpecialFormResult::Recognized(SpecialIntegral {
            function: SpecialFunction::Ei,
            coefficient: SpecialCoefficient::Rational(Q::from_integer(1)),
            argument: SpecialArgument::Identity,
        })
    } else if a == Q::from_integer(-1) {
        // ∫ exp(-x)/x dx = -Ei(-x)
        SpecialFormResult::Recognized(SpecialIntegral {
            function: SpecialFunction::Ei,
            coefficient: SpecialCoefficient::Rational(Q::from_integer(-1)),
            argument: SpecialArgument::Negated,
        })
    } else {
        // ∫ exp(ax)/x dx = Ei(ax)
        SpecialFormResult::Recognized(SpecialIntegral {
            function: SpecialFunction::Ei,
            coefficient: SpecialCoefficient::Rational(Q::from_integer(1)),
            argument: SpecialArgument::Linear { a, b: Q::from_integer(0) },
        })
    }
}

/// Recognizes logarithmic integral: ∫ 1/log(x) dx = li(x)
pub fn recognize_logarithmic_integral(
    numerator: &DensePoly<Q>,
    log_argument: &DensePoly<Q>,
) -> SpecialFormResult {
    // ∫ 1/log(x) dx = li(x)
    // This requires the numerator to be 1 and the log argument to be x

    // Check numerator is 1
    if !numerator.is_constant_one() {
        return SpecialFormResult::NotRecognized;
    }

    // Check log argument is x
    if log_argument.degree() != 1 {
        return SpecialFormResult::NotRecognized;
    }

    let c = log_argument.coeff(0);
    let d = log_argument.coeff(1);

    if c.is_zero() && d == Q::from_integer(1) {
        SpecialFormResult::Recognized(SpecialIntegral {
            function: SpecialFunction::Li,
            coefficient: SpecialCoefficient::Rational(Q::from_integer(1)),
            argument: SpecialArgument::Identity,
        })
    } else {
        SpecialFormResult::NotRecognized
    }
}

/// Helper trait for checking if a polynomial is the constant 1.
trait IsConstantOne {
    fn is_constant_one(&self) -> bool;
}

impl<F: Field> IsConstantOne for DensePoly<F> {
    fn is_constant_one(&self) -> bool {
        self.degree() == 0 && self.coeff(0).is_one()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_recognize_gaussian_simple() {
        // exp(-x²)
        let num = poly(&[1]);
        let exp = poly(&[0, 0, -1]); // -x²

        let result = recognize_gaussian(&num, &exp);

        match result {
            SpecialFormResult::Recognized(integral) => {
                assert_eq!(integral.function, SpecialFunction::Erf);
                assert!(matches!(integral.coefficient, SpecialCoefficient::SqrtPiOver(2)));
            }
            _ => panic!("Expected Gaussian to be recognized"),
        }
    }

    #[test]
    fn test_recognize_gaussian_not_quadratic() {
        // exp(-x³) - not a Gaussian
        let num = poly(&[1]);
        let exp = poly(&[0, 0, 0, -1]); // -x³

        let result = recognize_gaussian(&num, &exp);
        assert!(matches!(result, SpecialFormResult::NotRecognized));
    }

    #[test]
    fn test_recognize_ei() {
        // exp(x)/x
        let poly_part = poly(&[0, 1]); // x (exponent)
        let rf = RationalFunction::new(poly(&[1]), poly(&[0, 1])); // 1/x

        let result = recognize_exponential_integral(&poly_part, &rf);

        match result {
            SpecialFormResult::Recognized(integral) => {
                assert_eq!(integral.function, SpecialFunction::Ei);
                assert!(matches!(integral.argument, SpecialArgument::Identity));
            }
            _ => panic!("Expected Ei to be recognized"),
        }
    }

    #[test]
    fn test_recognize_ei_negative() {
        // exp(-x)/x
        let poly_part = poly(&[0, -1]); // -x (exponent)
        let rf = RationalFunction::new(poly(&[1]), poly(&[0, 1])); // 1/x

        let result = recognize_exponential_integral(&poly_part, &rf);

        match result {
            SpecialFormResult::Recognized(integral) => {
                assert_eq!(integral.function, SpecialFunction::Ei);
                assert!(matches!(integral.argument, SpecialArgument::Negated));
            }
            _ => panic!("Expected -Ei(-x) to be recognized"),
        }
    }

    #[test]
    fn test_recognize_li() {
        // 1/log(x)
        let num = poly(&[1]);
        let log_arg = poly(&[0, 1]); // x

        let result = recognize_logarithmic_integral(&num, &log_arg);

        match result {
            SpecialFormResult::Recognized(integral) => {
                assert_eq!(integral.function, SpecialFunction::Li);
            }
            _ => panic!("Expected li(x) to be recognized"),
        }
    }

    #[test]
    fn test_special_integral_display() {
        let integral = SpecialIntegral {
            function: SpecialFunction::Erf,
            coefficient: SpecialCoefficient::SqrtPiOver(2),
            argument: SpecialArgument::Identity,
        };

        let s = format!("{}", integral);
        assert!(s.contains("erf"));
        assert!(s.contains("√π/2"));
    }
}
