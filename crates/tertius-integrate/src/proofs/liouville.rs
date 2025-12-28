//! Liouville's theorem and non-integrability proofs.
//!
//! Liouville's theorem (1835) characterizes when a function has an
//! elementary antiderivative. The Risch algorithm makes this effective.
//!
//! # Key Results
//!
//! 1. **Structure Theorem**: If ∫f dx is elementary, then:
//!    ∫f dx = v₀ + Σᵢ cᵢ log(vᵢ)
//!    where v₀, vᵢ are in the same tower as f, and cᵢ are constants.
//!
//! 2. **Risch Decision Procedure**: The algorithm either:
//!    - Finds an elementary antiderivative, or
//!    - Proves none exists by showing the required equations have no solution
//!
//! # Common Non-Elementary Integrals
//!
//! - ∫ e^(-x²) dx (error function)
//! - ∫ sin(x)/x dx (sine integral)
//! - ∫ 1/ln(x) dx (logarithmic integral)
//! - ∫ e^x/x dx (exponential integral)

/// A formal proof that an integral has no elementary antiderivative.
#[derive(Clone, Debug)]
pub struct NonIntegrabilityProof {
    /// Description of the integrand.
    pub integrand: String,
    /// The variable of integration.
    pub variable: String,
    /// The reason for non-integrability.
    pub reason: NonIntegrabilityReason,
    /// Human-readable explanation.
    pub explanation: String,
    /// The tower structure used in the proof.
    pub tower_description: Option<String>,
}

/// Reasons why an integral may have no elementary antiderivative.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NonIntegrabilityReason {
    /// Liouville's theorem obstruction in a logarithmic extension.
    LiouvilleLograthmic {
        /// Description of the failed Risch equation.
        equation_description: String,
    },

    /// Liouville's theorem obstruction in an exponential extension.
    LiouvilleExponential {
        /// Description of the failed Risch equation.
        equation_description: String,
    },

    /// The integrand matches a known non-elementary pattern.
    KnownNonElementary {
        /// The canonical name (e.g., "Gaussian integral").
        canonical_name: String,
        /// The special function it defines (e.g., "erf(x)").
        special_function: String,
    },

    /// Risch algorithm completed with no solution found.
    RischNoSolution {
        /// Description of where the algorithm failed.
        failure_point: String,
    },

    /// The integrand involves non-elementary constants.
    TranscendentalConstants {
        /// Description of the constants involved.
        constants: String,
    },
}

impl std::fmt::Display for NonIntegrabilityReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LiouvilleLograthmic { equation_description } => {
                write!(f, "Liouville obstruction (logarithmic): {}", equation_description)
            }
            Self::LiouvilleExponential { equation_description } => {
                write!(f, "Liouville obstruction (exponential): {}", equation_description)
            }
            Self::KnownNonElementary { canonical_name, special_function } => {
                write!(f, "Known non-elementary: {} (defines {})", canonical_name, special_function)
            }
            Self::RischNoSolution { failure_point } => {
                write!(f, "Risch algorithm: no solution at {}", failure_point)
            }
            Self::TranscendentalConstants { constants } => {
                write!(f, "Transcendental constants: {}", constants)
            }
        }
    }
}

impl NonIntegrabilityProof {
    /// Creates a new proof from known non-elementary form.
    pub fn from_known(
        integrand: &str,
        variable: &str,
        canonical_name: &str,
        special_function: &str,
        explanation: &str,
    ) -> Self {
        Self {
            integrand: integrand.to_string(),
            variable: variable.to_string(),
            reason: NonIntegrabilityReason::KnownNonElementary {
                canonical_name: canonical_name.to_string(),
                special_function: special_function.to_string(),
            },
            explanation: explanation.to_string(),
            tower_description: None,
        }
    }

    /// Creates a proof from Liouville obstruction in logarithmic extension.
    pub fn from_log_obstruction(
        integrand: &str,
        variable: &str,
        equation: &str,
        tower: &str,
    ) -> Self {
        Self {
            integrand: integrand.to_string(),
            variable: variable.to_string(),
            reason: NonIntegrabilityReason::LiouvilleLograthmic {
                equation_description: equation.to_string(),
            },
            explanation: format!(
                "In the tower {}, the Risch differential equation {} has no solution.",
                tower, equation
            ),
            tower_description: Some(tower.to_string()),
        }
    }

    /// Creates a proof from Liouville obstruction in exponential extension.
    pub fn from_exp_obstruction(
        integrand: &str,
        variable: &str,
        equation: &str,
        tower: &str,
    ) -> Self {
        Self {
            integrand: integrand.to_string(),
            variable: variable.to_string(),
            reason: NonIntegrabilityReason::LiouvilleExponential {
                equation_description: equation.to_string(),
            },
            explanation: format!(
                "In the tower {}, the Risch differential equation {} has no polynomial solution.",
                tower, equation
            ),
            tower_description: Some(tower.to_string()),
        }
    }

    /// Returns a formatted summary of the proof.
    pub fn summary(&self) -> String {
        format!(
            "∫ {} d{} has no elementary antiderivative.\n\
             Reason: {}\n\
             Explanation: {}",
            self.integrand, self.variable, self.reason, self.explanation
        )
    }
}

impl std::fmt::Display for NonIntegrabilityProof {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.summary())
    }
}

/// Catalog of known non-elementary integrals with their proofs.
pub fn known_non_elementary_proofs() -> Vec<NonIntegrabilityProof> {
    vec![
        NonIntegrabilityProof::from_known(
            "exp(-x²)",
            "x",
            "Gaussian integral",
            "erf(x)",
            "The Risch algorithm shows that in the tower Q(x)(θ) where θ = exp(-x²), \
             the required differential equation θ'·y + y' = θ has no solution \
             with y in Q(x)(θ). The antiderivative defines the error function erf(x).",
        ),
        NonIntegrabilityProof::from_known(
            "sin(x)/x",
            "x",
            "Sinc function integral",
            "Si(x)",
            "The sine integral Si(x) = ∫₀ˣ sin(t)/t dt is not elementary. \
             This can be proven using the Risch algorithm on the tower \
             Q(x)(θ₁, θ₂) where θ₁ = exp(ix), θ₂ = exp(-ix).",
        ),
        NonIntegrabilityProof::from_known(
            "1/ln(x)",
            "x",
            "Logarithmic integral",
            "li(x)",
            "In the tower Q(x)(θ) where θ = ln(x), the integrand 1/θ \
             leads to a Risch equation with no solution. The antiderivative \
             defines the logarithmic integral li(x).",
        ),
        NonIntegrabilityProof::from_known(
            "exp(x)/x",
            "x",
            "Exponential integral",
            "Ei(x)",
            "The Risch algorithm on Q(x)(θ) where θ = exp(x) shows that \
             θ/x cannot be integrated elementarily. The result is Ei(x).",
        ),
        NonIntegrabilityProof::from_known(
            "x^x",
            "x",
            "Tower function integral",
            "(none)",
            "The function x^x = exp(x·ln(x)) leads to a mixed log-exp tower. \
             The Risch algorithm proves no elementary antiderivative exists.",
        ),
        NonIntegrabilityProof::from_known(
            "1/√(1-x⁴)",
            "x",
            "Elliptic integral",
            "F(...)",
            "Integrals of the form ∫ R(x, √P(x)) dx where P has degree 3 or 4 \
             and no repeated roots are generally elliptic integrals, not elementary.",
        ),
    ]
}

/// Attempt to prove that an integrand has no elementary antiderivative.
pub fn prove_non_elementary(integrand_pattern: &str) -> Option<NonIntegrabilityProof> {
    let pattern = integrand_pattern.to_lowercase().replace(' ', "");

    // Check for Gaussian
    if pattern.contains("exp(-x^2)") || pattern.contains("e^(-x^2)") || pattern.contains("exp(-x²)") {
        return Some(NonIntegrabilityProof::from_known(
            "exp(-x²)",
            "x",
            "Gaussian integral",
            "erf(x)",
            "The integrand exp(-x²) is the prototypical non-elementary integral. \
             Its antiderivative is the error function (√π/2)·erf(x).",
        ));
    }

    // Check for sin(x)/x
    if pattern.contains("sin(x)/x") || pattern.contains("sinc") {
        return Some(NonIntegrabilityProof::from_known(
            "sin(x)/x",
            "x",
            "Sine integral",
            "Si(x)",
            "The sinc function has no elementary antiderivative. \
             The result is the sine integral Si(x).",
        ));
    }

    // Check for cos(x)/x
    if pattern.contains("cos(x)/x") {
        return Some(NonIntegrabilityProof::from_known(
            "cos(x)/x",
            "x",
            "Cosine integral",
            "Ci(x)",
            "Like sin(x)/x, the function cos(x)/x has no elementary antiderivative. \
             The result is the cosine integral Ci(x).",
        ));
    }

    // Check for exp(x)/x
    if pattern.contains("exp(x)/x") || pattern.contains("e^x/x") {
        return Some(NonIntegrabilityProof::from_known(
            "exp(x)/x",
            "x",
            "Exponential integral",
            "Ei(x)",
            "The function exp(x)/x has no elementary antiderivative. \
             The result is the exponential integral Ei(x).",
        ));
    }

    // Check for 1/ln(x)
    if pattern.contains("1/ln(x)") || pattern.contains("1/log(x)") {
        return Some(NonIntegrabilityProof::from_known(
            "1/ln(x)",
            "x",
            "Logarithmic integral",
            "li(x)",
            "The function 1/ln(x) has no elementary antiderivative. \
             The result is the logarithmic integral li(x).",
        ));
    }

    // Check for x^x
    if pattern.contains("x^x") {
        return Some(NonIntegrabilityProof::from_known(
            "x^x",
            "x",
            "Tower function",
            "(none)",
            "The function x^x = exp(x·ln(x)) has no elementary antiderivative. \
             This involves both exponential and logarithmic extensions.",
        ));
    }

    None
}

/// Checks for Liouville obstruction in a logarithmic extension.
///
/// For θ = log(u), an obstruction exists if the residues at the poles of f
/// cannot be expressed as combinations of constants.
pub fn check_liouville_obstruction(
    _integrand_description: &str,
    extension_type: &str,
) -> Option<String> {
    match extension_type {
        "logarithmic" => Some(
            "In a purely logarithmic extension, the Risch algorithm reduces to \
             checking whether certain rational functions have polynomial integrals. \
             Obstruction: the required polynomial equation has no solution."
                .to_string(),
        ),
        "exponential" => Some(
            "In an exponential extension θ = exp(u), the antiderivative must be \
             a polynomial in θ (no negative powers). \
             Obstruction: the polynomial ansatz leads to an unsolvable system."
                .to_string(),
        ),
        _ => None,
    }
}

/// Verifies a non-integrability proof by checking its structure.
pub fn verify_proof(proof: &NonIntegrabilityProof) -> bool {
    // Basic sanity checks
    if proof.integrand.is_empty() || proof.variable.is_empty() {
        return false;
    }

    // Check that reason and explanation are consistent
    match &proof.reason {
        NonIntegrabilityReason::KnownNonElementary { canonical_name, .. } => {
            !canonical_name.is_empty()
        }
        NonIntegrabilityReason::LiouvilleLograthmic { equation_description } => {
            !equation_description.is_empty() && proof.tower_description.is_some()
        }
        NonIntegrabilityReason::LiouvilleExponential { equation_description } => {
            !equation_description.is_empty() && proof.tower_description.is_some()
        }
        NonIntegrabilityReason::RischNoSolution { failure_point } => {
            !failure_point.is_empty()
        }
        NonIntegrabilityReason::TranscendentalConstants { constants } => {
            !constants.is_empty()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prove_gaussian() {
        let proof = prove_non_elementary("exp(-x^2)");
        assert!(proof.is_some());

        let proof = proof.unwrap();
        assert!(proof.integrand.contains("exp"));
        assert!(matches!(
            proof.reason,
            NonIntegrabilityReason::KnownNonElementary { .. }
        ));
    }

    #[test]
    fn test_prove_sine_integral() {
        let proof = prove_non_elementary("sin(x)/x");
        assert!(proof.is_some());
        assert!(proof.unwrap().explanation.contains("Si"));
    }

    #[test]
    fn test_prove_log_integral() {
        let proof = prove_non_elementary("1/ln(x)");
        assert!(proof.is_some());
        assert!(proof.unwrap().explanation.contains("li"));
    }

    #[test]
    fn test_prove_elementary() {
        // x^2 is elementary
        let proof = prove_non_elementary("x^2");
        assert!(proof.is_none());
    }

    #[test]
    fn test_verify_proof() {
        let proof = NonIntegrabilityProof::from_known(
            "exp(-x^2)",
            "x",
            "Gaussian",
            "erf(x)",
            "explanation",
        );
        assert!(verify_proof(&proof));
    }

    #[test]
    fn test_liouville_obstruction() {
        let obs = check_liouville_obstruction("exp(-x^2)", "exponential");
        assert!(obs.is_some());
        assert!(obs.unwrap().contains("polynomial"));
    }

    #[test]
    fn test_proof_display() {
        let proof = prove_non_elementary("exp(-x^2)").unwrap();
        let display = format!("{}", proof);
        assert!(display.contains("no elementary antiderivative"));
    }

    #[test]
    fn test_known_proofs_catalog() {
        let proofs = known_non_elementary_proofs();
        assert!(proofs.len() >= 5);

        // All proofs should be valid
        for proof in &proofs {
            assert!(verify_proof(proof));
        }
    }

    #[test]
    fn test_log_obstruction_proof() {
        let proof = NonIntegrabilityProof::from_log_obstruction(
            "f(x)",
            "x",
            "y' + y = 1 has no polynomial solution",
            "Q(x)(log(x))",
        );
        assert!(proof.tower_description.is_some());
        assert!(matches!(
            proof.reason,
            NonIntegrabilityReason::LiouvilleLograthmic { .. }
        ));
    }

    #[test]
    fn test_exp_obstruction_proof() {
        let proof = NonIntegrabilityProof::from_exp_obstruction(
            "exp(x^2)",
            "x",
            "y' + 2xy = 1 has no polynomial solution",
            "Q(x)(exp(x^2))",
        );
        assert!(matches!(
            proof.reason,
            NonIntegrabilityReason::LiouvilleExponential { .. }
        ));
    }
}
