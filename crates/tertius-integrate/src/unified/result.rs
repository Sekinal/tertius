//! Integration result types.
//!
//! Defines the unified `IntegrationResult` enum that captures all possible
//! outcomes from various integration methods.

use tertius_core::handle::ExprHandle;
use tertius_rings::rationals::Q;

/// Unified result type for all integration methods.
#[derive(Clone, Debug)]
pub enum IntegrationResult {
    /// A closed-form symbolic antiderivative was found.
    Symbolic(SymbolicAntiderivative),

    /// The integral involves special functions (erf, Ei, Li, etc.).
    SpecialFunction(SpecialFunctionResult),

    /// The integral is expressible in terms of elliptic integrals.
    Elliptic(EllipticResult),

    /// A definite integral was computed numerically.
    Numerical(NumericalResult),

    /// No elementary antiderivative exists (with proof).
    NonElementary(NonElementaryResult),

    /// The integration method could not determine a result.
    Unknown(UnknownReason),
}

impl IntegrationResult {
    /// Returns true if a symbolic result was found.
    pub fn is_symbolic(&self) -> bool {
        matches!(self, IntegrationResult::Symbolic(_))
    }

    /// Returns true if the result involves special functions.
    pub fn is_special_function(&self) -> bool {
        matches!(self, IntegrationResult::SpecialFunction(_))
    }

    /// Returns true if the integral is non-elementary.
    pub fn is_non_elementary(&self) -> bool {
        matches!(self, IntegrationResult::NonElementary(_))
    }

    /// Returns true if a result was found (symbolic, special, elliptic, or numerical).
    pub fn is_success(&self) -> bool {
        !matches!(self, IntegrationResult::Unknown(_))
    }

    /// Extracts the symbolic antiderivative if available.
    pub fn as_symbolic(&self) -> Option<&SymbolicAntiderivative> {
        match self {
            IntegrationResult::Symbolic(s) => Some(s),
            _ => None,
        }
    }

    /// Extracts the numerical result if available.
    pub fn as_numerical(&self) -> Option<&NumericalResult> {
        match self {
            IntegrationResult::Numerical(n) => Some(n),
            _ => None,
        }
    }
}

/// Symbolic antiderivative in closed form.
#[derive(Clone, Debug)]
pub struct SymbolicAntiderivative {
    /// The antiderivative expression.
    pub result: ExprHandle,

    /// Integration method used.
    pub method: IntegrationMethod,

    /// Whether verification (differentiation) succeeded.
    pub verified: bool,

    /// String representation for display.
    pub display: String,
}

impl SymbolicAntiderivative {
    /// Creates a new symbolic antiderivative.
    pub fn new(result: ExprHandle, method: IntegrationMethod) -> Self {
        Self {
            result,
            method,
            verified: false,
            display: String::new(),
        }
    }

    /// Sets the verified flag.
    pub fn with_verified(mut self, verified: bool) -> Self {
        self.verified = verified;
        self
    }

    /// Sets the display string.
    pub fn with_display(mut self, display: String) -> Self {
        self.display = display;
        self
    }
}

/// Method used for integration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntegrationMethod {
    /// Simple power rule: ∫x^n dx = x^(n+1)/(n+1)
    PowerRule,
    /// Hermite reduction + Rothstein-Trager for rational functions
    RationalHermiteRT,
    /// Risch algorithm for logarithmic extensions
    RischLogarithmic,
    /// Risch algorithm for exponential extensions
    RischExponential,
    /// Rationalization of algebraic functions
    Rationalization,
    /// Table lookup / pattern matching
    TableLookup,
    /// Heuristic methods
    Heuristic,
    /// Integration by parts
    ByParts,
    /// Substitution method
    Substitution,
}

impl std::fmt::Display for IntegrationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IntegrationMethod::PowerRule => write!(f, "power rule"),
            IntegrationMethod::RationalHermiteRT => write!(f, "Hermite + Rothstein-Trager"),
            IntegrationMethod::RischLogarithmic => write!(f, "Risch (logarithmic)"),
            IntegrationMethod::RischExponential => write!(f, "Risch (exponential)"),
            IntegrationMethod::Rationalization => write!(f, "rationalization"),
            IntegrationMethod::TableLookup => write!(f, "table lookup"),
            IntegrationMethod::Heuristic => write!(f, "heuristic"),
            IntegrationMethod::ByParts => write!(f, "integration by parts"),
            IntegrationMethod::Substitution => write!(f, "substitution"),
        }
    }
}

/// Result involving special functions.
#[derive(Clone, Debug)]
pub struct SpecialFunctionResult {
    /// The special function(s) involved.
    pub functions: Vec<SpecialFunctionTerm>,

    /// String representation for display.
    pub display: String,

    /// The result expression handle (if available in arena).
    pub result: Option<ExprHandle>,
}

impl SpecialFunctionResult {
    /// Creates a new special function result.
    pub fn new(display: String) -> Self {
        Self {
            functions: Vec::new(),
            display,
            result: None,
        }
    }

    /// Adds a special function term.
    pub fn with_term(mut self, term: SpecialFunctionTerm) -> Self {
        self.functions.push(term);
        self
    }
}

/// A term involving a special function.
#[derive(Clone, Debug)]
pub struct SpecialFunctionTerm {
    /// The special function.
    pub function: SpecialFunction,
    /// Coefficient (rational).
    pub coefficient: Q,
    /// Argument description.
    pub argument: String,
}

/// Special functions that appear in integration results.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpecialFunction {
    /// Error function erf(x)
    Erf,
    /// Complementary error function erfc(x)
    Erfc,
    /// Exponential integral Ei(x)
    Ei,
    /// Logarithmic integral li(x)
    Li,
    /// Sine integral Si(x)
    Si,
    /// Cosine integral Ci(x)
    Ci,
    /// Fresnel sine integral S(x)
    FresnelS,
    /// Fresnel cosine integral C(x)
    FresnelC,
    /// Gamma function Γ(x)
    Gamma,
    /// Incomplete gamma function γ(s, x)
    IncompleteGamma,
    /// Dilogarithm Li₂(x)
    Dilog,
    /// Polylogarithm Li_n(x)
    Polylog,
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
            SpecialFunction::IncompleteGamma => write!(f, "γ"),
            SpecialFunction::Dilog => write!(f, "Li₂"),
            SpecialFunction::Polylog => write!(f, "Li"),
        }
    }
}

/// Elliptic integral result.
#[derive(Clone, Debug)]
pub struct EllipticResult {
    /// The kind of elliptic integral.
    pub kind: EllipticKind,
    /// Amplitude (as expression handle).
    pub amplitude: Option<ExprHandle>,
    /// Modulus k.
    pub modulus: f64,
    /// Parameter n (for third kind).
    pub parameter: Option<f64>,
    /// String representation.
    pub display: String,
}

impl EllipticResult {
    /// Creates a new elliptic result.
    pub fn new(kind: EllipticKind, modulus: f64) -> Self {
        Self {
            kind,
            amplitude: None,
            modulus,
            parameter: None,
            display: String::new(),
        }
    }
}

/// Kinds of elliptic integrals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EllipticKind {
    /// Elliptic integral of the first kind F(φ, k)
    FirstKind,
    /// Elliptic integral of the second kind E(φ, k)
    SecondKind,
    /// Elliptic integral of the third kind Π(n; φ, k)
    ThirdKind,
}

impl std::fmt::Display for EllipticKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EllipticKind::FirstKind => write!(f, "F"),
            EllipticKind::SecondKind => write!(f, "E"),
            EllipticKind::ThirdKind => write!(f, "Π"),
        }
    }
}

/// Numerical integration result.
#[derive(Clone, Debug)]
pub struct NumericalResult {
    /// The computed value.
    pub value: f64,
    /// Error estimate.
    pub error_estimate: f64,
    /// Number of function evaluations.
    pub evaluations: usize,
    /// Whether convergence was achieved.
    pub converged: bool,
}

impl NumericalResult {
    /// Creates a new numerical result.
    pub fn new(value: f64, error_estimate: f64) -> Self {
        Self {
            value,
            error_estimate,
            evaluations: 0,
            converged: true,
        }
    }

    /// Sets the number of evaluations.
    pub fn with_evaluations(mut self, evaluations: usize) -> Self {
        self.evaluations = evaluations;
        self
    }

    /// Sets the convergence flag.
    pub fn with_converged(mut self, converged: bool) -> Self {
        self.converged = converged;
        self
    }
}

impl std::fmt::Display for NumericalResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.converged {
            write!(f, "{:.10} ± {:.2e}", self.value, self.error_estimate)
        } else {
            write!(
                f,
                "{:.10} ± {:.2e} (not converged)",
                self.value, self.error_estimate
            )
        }
    }
}

/// Non-elementary integral with proof.
#[derive(Clone, Debug)]
pub struct NonElementaryResult {
    /// Reason for non-integrability.
    pub reason: NonElementaryReason,
    /// Proof or explanation.
    pub proof: Option<String>,
    /// Canonical special function form if known.
    pub canonical_form: Option<String>,
}

impl NonElementaryResult {
    /// Creates a new non-elementary result.
    pub fn new(reason: NonElementaryReason) -> Self {
        Self {
            reason,
            proof: None,
            canonical_form: None,
        }
    }

    /// Sets the proof.
    pub fn with_proof(mut self, proof: String) -> Self {
        self.proof = Some(proof);
        self
    }

    /// Sets the canonical form.
    pub fn with_canonical(mut self, form: String) -> Self {
        self.canonical_form = Some(form);
        self
    }
}

/// Reason why an integral is non-elementary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NonElementaryReason {
    /// Liouville's theorem (logarithmic case)
    LiouvilleLogarithmic,
    /// Liouville's theorem (exponential case)
    LiouvilleExponential,
    /// Risch algorithm decision procedure
    RischDecision,
    /// Known non-elementary form
    KnownNonElementary(String),
}

impl std::fmt::Display for NonElementaryReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NonElementaryReason::LiouvilleLogarithmic => {
                write!(f, "Liouville's theorem (logarithmic obstruction)")
            }
            NonElementaryReason::LiouvilleExponential => {
                write!(f, "Liouville's theorem (exponential obstruction)")
            }
            NonElementaryReason::RischDecision => write!(f, "Risch decision procedure"),
            NonElementaryReason::KnownNonElementary(name) => {
                write!(f, "known non-elementary: {}", name)
            }
        }
    }
}

/// Reason why integration could not determine a result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum UnknownReason {
    /// Computation timed out
    Timeout,
    /// Expression too complex
    ComplexityLimit,
    /// Unsupported expression form
    UnsupportedForm,
    /// Internal error
    InternalError(String),
}

impl std::fmt::Display for UnknownReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            UnknownReason::Timeout => write!(f, "computation timed out"),
            UnknownReason::ComplexityLimit => write!(f, "expression too complex"),
            UnknownReason::UnsupportedForm => write!(f, "unsupported expression form"),
            UnknownReason::InternalError(msg) => write!(f, "internal error: {}", msg),
        }
    }
}

impl std::fmt::Display for IntegrationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IntegrationResult::Symbolic(s) => {
                if s.display.is_empty() {
                    write!(f, "symbolic result ({})", s.method)
                } else {
                    write!(f, "{}", s.display)
                }
            }
            IntegrationResult::SpecialFunction(sf) => write!(f, "{}", sf.display),
            IntegrationResult::Elliptic(e) => write!(f, "{}", e.display),
            IntegrationResult::Numerical(n) => write!(f, "{}", n),
            IntegrationResult::NonElementary(ne) => {
                if let Some(ref form) = ne.canonical_form {
                    write!(f, "non-elementary: {} (canonical: {})", ne.reason, form)
                } else {
                    write!(f, "non-elementary: {}", ne.reason)
                }
            }
            IntegrationResult::Unknown(reason) => write!(f, "unknown: {}", reason),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integration_method_display() {
        assert_eq!(format!("{}", IntegrationMethod::PowerRule), "power rule");
        assert_eq!(
            format!("{}", IntegrationMethod::RationalHermiteRT),
            "Hermite + Rothstein-Trager"
        );
    }

    #[test]
    fn test_numerical_result_display() {
        let result = NumericalResult::new(3.14159, 1e-10);
        let display = format!("{}", result);
        assert!(display.contains("3.14159"));
    }

    #[test]
    fn test_unknown_reason_display() {
        assert_eq!(format!("{}", UnknownReason::Timeout), "computation timed out");
        assert_eq!(
            format!("{}", UnknownReason::InternalError("test".to_string())),
            "internal error: test"
        );
    }
}
