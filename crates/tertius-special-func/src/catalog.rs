//! Catalog of special functions and pattern recognition.
//!
//! This module provides:
//! - A unified representation of special functions
//! - Pattern matching to recognize integrands that reduce to special functions
//! - Lookup of known non-elementary integrals

use crate::error_func::ErrorFunction;
use crate::polylog::Polylogarithm;
use crate::elliptic::{EllipticE, EllipticF};
use crate::hypergeometric::Hypergeometric2F1;

/// A special function that may arise from integration.
#[derive(Clone, Debug)]
pub enum SpecialFunction {
    /// Polylogarithm Li_n(x)
    Polylog(Polylogarithm),
    /// Error function erf(x) or variants
    ErrorFunc(ErrorFunction),
    /// Elliptic integral of first kind F(φ, k)
    EllipticF(EllipticF),
    /// Elliptic integral of second kind E(φ, k)
    EllipticE(EllipticE),
    /// Hypergeometric ₂F₁(a, b; c; z)
    Hypergeometric(Hypergeometric2F1),
    /// Logarithmic integral li(x) = ∫₀ˣ dt/ln(t)
    LogarithmicIntegral { argument: String },
    /// Sine integral Si(x) = ∫₀ˣ sin(t)/t dt
    SineIntegral { argument: String },
    /// Cosine integral Ci(x) = -∫ₓ^∞ cos(t)/t dt
    CosineIntegral { argument: String },
    /// Exponential integral Ei(x) = ∫₋∞ˣ eᵗ/t dt
    ExponentialIntegral { argument: String },
    /// Fresnel sine integral S(x) = ∫₀ˣ sin(πt²/2) dt
    FresnelS { argument: String },
    /// Fresnel cosine integral C(x) = ∫₀ˣ cos(πt²/2) dt
    FresnelC { argument: String },
}

impl std::fmt::Display for SpecialFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpecialFunction::Polylog(p) => write!(f, "{}", p),
            SpecialFunction::ErrorFunc(e) => write!(f, "{}", e),
            SpecialFunction::EllipticF(ef) => write!(f, "{}", ef),
            SpecialFunction::EllipticE(ee) => write!(f, "{}", ee),
            SpecialFunction::Hypergeometric(h) => write!(f, "{}", h),
            SpecialFunction::LogarithmicIntegral { argument } => write!(f, "li({})", argument),
            SpecialFunction::SineIntegral { argument } => write!(f, "Si({})", argument),
            SpecialFunction::CosineIntegral { argument } => write!(f, "Ci({})", argument),
            SpecialFunction::ExponentialIntegral { argument } => write!(f, "Ei({})", argument),
            SpecialFunction::FresnelS { argument } => write!(f, "S({})", argument),
            SpecialFunction::FresnelC { argument } => write!(f, "C({})", argument),
        }
    }
}

/// Known non-elementary integrals with their special function results.
#[derive(Clone, Debug)]
pub struct NonElementaryIntegral {
    /// Description of the integrand.
    pub integrand: String,
    /// The special function result.
    pub result: SpecialFunction,
    /// Human-readable explanation.
    pub explanation: String,
}

/// Catalog of known non-elementary integrals.
pub fn known_non_elementary_integrals() -> Vec<NonElementaryIntegral> {
    vec![
        NonElementaryIntegral {
            integrand: "exp(-x²)".to_string(),
            result: SpecialFunction::ErrorFunc(ErrorFunction::erf("x")),
            explanation: "∫ exp(-x²) dx = (√π/2) erf(x) + C".to_string(),
        },
        NonElementaryIntegral {
            integrand: "1/ln(x)".to_string(),
            result: SpecialFunction::LogarithmicIntegral {
                argument: "x".to_string(),
            },
            explanation: "∫ 1/ln(x) dx = li(x) + C".to_string(),
        },
        NonElementaryIntegral {
            integrand: "sin(x)/x".to_string(),
            result: SpecialFunction::SineIntegral {
                argument: "x".to_string(),
            },
            explanation: "∫ sin(x)/x dx = Si(x) + C".to_string(),
        },
        NonElementaryIntegral {
            integrand: "cos(x)/x".to_string(),
            result: SpecialFunction::CosineIntegral {
                argument: "x".to_string(),
            },
            explanation: "∫ cos(x)/x dx = Ci(x) + C".to_string(),
        },
        NonElementaryIntegral {
            integrand: "exp(x)/x".to_string(),
            result: SpecialFunction::ExponentialIntegral {
                argument: "x".to_string(),
            },
            explanation: "∫ exp(x)/x dx = Ei(x) + C".to_string(),
        },
        NonElementaryIntegral {
            integrand: "ln(x)/(x-1)".to_string(),
            result: SpecialFunction::Polylog(Polylogarithm::dilog("x")),
            explanation: "∫ ln(x)/(x-1) dx = Li₂(x) + C".to_string(),
        },
        NonElementaryIntegral {
            integrand: "sin(x²)".to_string(),
            result: SpecialFunction::FresnelS {
                argument: "x".to_string(),
            },
            explanation: "∫ sin(x²) dx is related to Fresnel S".to_string(),
        },
        NonElementaryIntegral {
            integrand: "cos(x²)".to_string(),
            result: SpecialFunction::FresnelC {
                argument: "x".to_string(),
            },
            explanation: "∫ cos(x²) dx is related to Fresnel C".to_string(),
        },
    ]
}

/// Attempt to recognize an integrand pattern and return the special function result.
pub fn recognize_special_function(integrand_pattern: &str) -> Option<SpecialFunction> {
    let pattern = integrand_pattern.to_lowercase().replace(' ', "");

    // Check for Gaussian/error function patterns
    if pattern.contains("exp(-x^2)") || pattern.contains("e^(-x^2)") {
        return Some(SpecialFunction::ErrorFunc(ErrorFunction::erf("x")));
    }

    if pattern.contains("exp(-x²)") || pattern.contains("e^(-x²)") {
        return Some(SpecialFunction::ErrorFunc(ErrorFunction::erf("x")));
    }

    // Check for logarithmic integral
    if pattern.contains("1/ln(x)") || pattern.contains("1/log(x)") {
        return Some(SpecialFunction::LogarithmicIntegral {
            argument: "x".to_string(),
        });
    }

    // Check for sine integral
    if pattern.contains("sin(x)/x") || pattern.contains("sinc(x)") {
        return Some(SpecialFunction::SineIntegral {
            argument: "x".to_string(),
        });
    }

    // Check for cosine integral
    if pattern.contains("cos(x)/x") {
        return Some(SpecialFunction::CosineIntegral {
            argument: "x".to_string(),
        });
    }

    // Check for exponential integral
    if pattern.contains("exp(x)/x") || pattern.contains("e^x/x") {
        return Some(SpecialFunction::ExponentialIntegral {
            argument: "x".to_string(),
        });
    }

    // Check for dilogarithm pattern
    if pattern.contains("ln(x)/(x-1)") || pattern.contains("log(x)/(x-1)") {
        return Some(SpecialFunction::Polylog(Polylogarithm::dilog("x")));
    }

    if pattern.contains("ln(1-x)/x") {
        return Some(SpecialFunction::Polylog(Polylogarithm::dilog("1-x")));
    }

    // Check for Fresnel patterns
    if pattern.contains("sin(x^2)") || pattern.contains("sin(x²)") {
        return Some(SpecialFunction::FresnelS {
            argument: "x".to_string(),
        });
    }

    if pattern.contains("cos(x^2)") || pattern.contains("cos(x²)") {
        return Some(SpecialFunction::FresnelC {
            argument: "x".to_string(),
        });
    }

    None
}

/// Checks if an integrand is known to be non-elementary.
pub fn is_known_non_elementary(integrand_pattern: &str) -> bool {
    recognize_special_function(integrand_pattern).is_some()
}

/// Returns a human-readable explanation for why an integral is non-elementary.
pub fn explain_non_elementary(integrand_pattern: &str) -> Option<String> {
    if let Some(sf) = recognize_special_function(integrand_pattern) {
        let explanation = match sf {
            SpecialFunction::ErrorFunc(_) => {
                "The integrand involves exp(-x²) which has no elementary antiderivative. \
                 The result is expressed using the error function erf(x)."
            }
            SpecialFunction::LogarithmicIntegral { .. } => {
                "The integrand 1/ln(x) has no elementary antiderivative. \
                 The result is the logarithmic integral li(x)."
            }
            SpecialFunction::SineIntegral { .. } => {
                "The integrand sin(x)/x has no elementary antiderivative. \
                 The result is the sine integral Si(x)."
            }
            SpecialFunction::CosineIntegral { .. } => {
                "The integrand cos(x)/x has no elementary antiderivative. \
                 The result is the cosine integral Ci(x)."
            }
            SpecialFunction::ExponentialIntegral { .. } => {
                "The integrand exp(x)/x has no elementary antiderivative. \
                 The result is the exponential integral Ei(x)."
            }
            SpecialFunction::Polylog(_) => {
                "This integral leads to a polylogarithm Li_n(x), which is a special function."
            }
            SpecialFunction::FresnelS { .. } | SpecialFunction::FresnelC { .. } => {
                "Integrals of sin(x²) and cos(x²) lead to Fresnel integrals S(x) and C(x)."
            }
            _ => {
                "This integral has no elementary antiderivative and requires special functions."
            }
        };
        Some(explanation.to_string())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_recognize_gaussian() {
        let result = recognize_special_function("exp(-x^2)");
        assert!(matches!(result, Some(SpecialFunction::ErrorFunc(_))));
    }

    #[test]
    fn test_recognize_log_integral() {
        let result = recognize_special_function("1/ln(x)");
        assert!(matches!(result, Some(SpecialFunction::LogarithmicIntegral { .. })));
    }

    #[test]
    fn test_recognize_sine_integral() {
        let result = recognize_special_function("sin(x)/x");
        assert!(matches!(result, Some(SpecialFunction::SineIntegral { .. })));
    }

    #[test]
    fn test_recognize_dilog() {
        let result = recognize_special_function("ln(x)/(x-1)");
        assert!(matches!(result, Some(SpecialFunction::Polylog(_))));
    }

    #[test]
    fn test_known_non_elementary() {
        assert!(is_known_non_elementary("exp(-x^2)"));
        assert!(is_known_non_elementary("sin(x)/x"));
        assert!(!is_known_non_elementary("x^2"));
    }

    #[test]
    fn test_explain_non_elementary() {
        let explanation = explain_non_elementary("exp(-x^2)");
        assert!(explanation.is_some());
        assert!(explanation.unwrap().contains("error function"));
    }

    #[test]
    fn test_catalog_size() {
        let catalog = known_non_elementary_integrals();
        assert!(catalog.len() >= 8);
    }

    #[test]
    fn test_special_function_display() {
        let sf = SpecialFunction::SineIntegral {
            argument: "x".to_string(),
        };
        assert_eq!(format!("{}", sf), "Si(x)");
    }
}
