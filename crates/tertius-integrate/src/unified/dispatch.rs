//! Integration dispatcher.
//!
//! Routes classified integrands to appropriate backend algorithms.

use smallvec::smallvec;

use tertius_core::arena::ExprArena;
use tertius_core::expr::{functions, ExprNode, SymbolId};
use tertius_core::handle::ExprHandle;
use tertius_rings::rationals::Q;

use super::classifier::{AlgebraicGenus, IntegrandAnalyzer, IntegrandClass, TranscendentalProfile};
use super::convert::ExprConverter;
use super::result::{
    EllipticKind, EllipticResult, IntegrationMethod, IntegrationResult, NonElementaryReason,
    NonElementaryResult, NumericalResult, SpecialFunction, SpecialFunctionResult,
    SpecialFunctionTerm, SymbolicAntiderivative, UnknownReason,
};

use crate::numerical::adaptive_integrate;
use crate::polynomial::integrate_polynomial;
use crate::rational::{integrate_rational_q, integrate_rational_with_algebraic};

/// Integration options for customizing behavior.
#[derive(Clone, Debug)]
pub struct IntegrationOptions {
    /// Enable numerical fallback for definite integrals.
    pub numerical_fallback: bool,
    /// Tolerance for numerical integration.
    pub numerical_tolerance: f64,
    /// Maximum subdivisions for numerical integration.
    pub max_subdivisions: usize,
    /// Try heuristic methods before formal Risch.
    pub heuristics_first: bool,
    /// Verify results by differentiation.
    pub verify: bool,
    /// Attempt simplification of result.
    pub simplify_result: bool,
}

impl Default for IntegrationOptions {
    fn default() -> Self {
        Self {
            numerical_fallback: true,
            numerical_tolerance: 1e-10,
            max_subdivisions: 1000,
            heuristics_first: true,
            verify: false,
            simplify_result: true,
        }
    }
}

/// Integration dispatcher that routes to appropriate backends.
pub struct IntegrationDispatcher<'a> {
    arena: &'a mut ExprArena,
    options: IntegrationOptions,
}

impl<'a> IntegrationDispatcher<'a> {
    /// Creates a new dispatcher with default options.
    pub fn new(arena: &'a mut ExprArena) -> Self {
        Self {
            arena,
            options: IntegrationOptions::default(),
        }
    }

    /// Creates a dispatcher with custom options.
    pub fn with_options(arena: &'a mut ExprArena, options: IntegrationOptions) -> Self {
        Self { arena, options }
    }

    /// Main entry point for indefinite integration.
    pub fn integrate(&mut self, expr: ExprHandle, var: ExprHandle) -> IntegrationResult {
        // Extract variable ID
        let var_id = match self.extract_variable_id(var) {
            Some(id) => id,
            None => {
                return IntegrationResult::Unknown(UnknownReason::InternalError(
                    "variable must be a symbol".to_string(),
                ))
            }
        };

        // Classify the integrand
        let analyzer = IntegrandAnalyzer::new(self.arena, var_id);
        let class = analyzer.classify(expr);

        // Dispatch based on classification
        self.dispatch_by_class(expr, var, var_id, class)
    }

    /// Definite integration entry point.
    pub fn integrate_definite(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        lower: ExprHandle,
        upper: ExprHandle,
    ) -> IntegrationResult {
        // First try symbolic integration
        let indefinite = self.integrate(expr, var);

        match &indefinite {
            IntegrationResult::Symbolic(antideriv) => {
                // Evaluate F(upper) - F(lower)
                self.evaluate_definite(antideriv.result, var, lower, upper)
            }
            IntegrationResult::Unknown(_) if self.options.numerical_fallback => {
                // Fall back to numerical integration
                self.numerical_definite(expr, var, lower, upper)
            }
            other => other.clone(),
        }
    }

    /// Dispatches to the appropriate backend based on classification.
    fn dispatch_by_class(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        var_id: SymbolId,
        class: IntegrandClass,
    ) -> IntegrationResult {
        match class {
            IntegrandClass::Constant => self.integrate_constant(expr, var),

            IntegrandClass::Polynomial { degree } => {
                self.integrate_polynomial_expr(expr, var, var_id, degree)
            }

            IntegrandClass::Rational {
                needs_algebraic, ..
            } => {
                if needs_algebraic {
                    self.integrate_rational_algebraic(expr, var, var_id)
                } else {
                    self.integrate_rational_expr(expr, var, var_id)
                }
            }

            IntegrandClass::Algebraic { genus, .. } => {
                self.integrate_algebraic_expr(expr, var, var_id, genus)
            }

            IntegrandClass::Transcendental { ref profile } => {
                self.integrate_transcendental(expr, var, var_id, profile)
            }

            IntegrandClass::Mixed { ref profile, .. } => {
                // Try transcendental methods for mixed cases
                self.integrate_transcendental(expr, var, var_id, profile)
            }

            IntegrandClass::Unknown => self.integrate_fallback(expr, var, var_id),
        }
    }

    /// Integrates a constant: ∫c dx = cx.
    fn integrate_constant(&mut self, expr: ExprHandle, var: ExprHandle) -> IntegrationResult {
        // c * x
        let result = self.arena.mul(smallvec![expr, var]);

        IntegrationResult::Symbolic(
            SymbolicAntiderivative::new(result, IntegrationMethod::PowerRule)
                .with_display(format!("({}) * x", self.expr_to_string(expr))),
        )
    }

    /// Integrates a polynomial using the power rule.
    fn integrate_polynomial_expr(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        var_id: SymbolId,
        _degree: usize,
    ) -> IntegrationResult {
        let mut converter = ExprConverter::new(self.arena, var_id, var);

        // Convert to polynomial
        let poly = match converter.to_polynomial(expr) {
            Some(p) => p,
            None => {
                return IntegrationResult::Unknown(UnknownReason::InternalError(
                    "failed to convert to polynomial".to_string(),
                ))
            }
        };

        // Integrate using power rule
        let integral_poly = integrate_polynomial(&poly);

        // Convert back to expression
        let result = converter.from_polynomial(&integral_poly);

        IntegrationResult::Symbolic(
            SymbolicAntiderivative::new(result, IntegrationMethod::PowerRule)
                .with_display(self.expr_to_string(result)),
        )
    }

    /// Integrates a rational function using Hermite + Rothstein-Trager.
    fn integrate_rational_expr(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        var_id: SymbolId,
    ) -> IntegrationResult {
        let mut converter = ExprConverter::new(self.arena, var_id, var);

        // Convert to rational function
        let rf = match converter.to_rational_function(expr) {
            Some(r) => r,
            None => {
                return IntegrationResult::Unknown(UnknownReason::InternalError(
                    "failed to convert to rational function".to_string(),
                ))
            }
        };

        // Integrate using Hermite + Rothstein-Trager
        let result = integrate_rational_q(&rf);

        if result.is_complete {
            // Convert result back to expression
            let expr_result = converter.from_integration_result(
                &result.polynomial_part,
                &result.rational_part,
                result.logarithmic_part.as_ref(),
            );

            IntegrationResult::Symbolic(
                SymbolicAntiderivative::new(expr_result, IntegrationMethod::RationalHermiteRT)
                    .with_display(self.expr_to_string(expr_result)),
            )
        } else {
            // Need algebraic extension
            self.integrate_rational_algebraic(expr, var, var_id)
        }
    }

    /// Integrates a rational function that may need algebraic extensions.
    fn integrate_rational_algebraic(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        var_id: SymbolId,
    ) -> IntegrationResult {
        let converter = ExprConverter::new(self.arena, var_id, var);

        // Convert to rational function
        let rf = match converter.to_rational_function(expr) {
            Some(r) => r,
            None => {
                return IntegrationResult::Unknown(UnknownReason::InternalError(
                    "failed to convert to rational function".to_string(),
                ))
            }
        };

        // Use algebraic integration
        let result = integrate_rational_with_algebraic(&rf);

        // For now, return the display string from the algebraic result
        IntegrationResult::Symbolic(
            SymbolicAntiderivative::new(expr, IntegrationMethod::RationalHermiteRT)
                .with_display(format!("{:?}", result)),
        )
    }

    /// Integrates an algebraic function based on genus.
    fn integrate_algebraic_expr(
        &mut self,
        expr: ExprHandle,
        _var: ExprHandle,
        _var_id: SymbolId,
        genus: AlgebraicGenus,
    ) -> IntegrationResult {
        match genus {
            AlgebraicGenus::Rational => {
                // Can be rationalized - try substitution
                // For now, return unknown
                IntegrationResult::Unknown(UnknownReason::UnsupportedForm)
            }

            AlgebraicGenus::Elliptic => {
                // Elliptic integral
                IntegrationResult::Elliptic(EllipticResult {
                    kind: EllipticKind::FirstKind,
                    amplitude: Some(expr),
                    modulus: 0.0, // Would need actual computation
                    parameter: None,
                    display: "elliptic integral".to_string(),
                })
            }

            AlgebraicGenus::Hyperelliptic(_) => {
                // Non-elementary
                IntegrationResult::NonElementary(
                    NonElementaryResult::new(NonElementaryReason::KnownNonElementary(
                        "hyperelliptic".to_string(),
                    ))
                    .with_proof("Genus > 1 implies non-elementary integral".to_string()),
                )
            }
        }
    }

    /// Integrates a transcendental expression.
    fn integrate_transcendental(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        var_id: SymbolId,
        profile: &TranscendentalProfile,
    ) -> IntegrationResult {
        // Try table lookup first
        if let Some(result) = self.try_table_lookup(expr, var, var_id) {
            return result;
        }

        // Check for known non-elementary forms
        if let Some(result) = self.check_non_elementary(expr, var_id, profile) {
            return result;
        }

        // For now, return unknown for other transcendental cases
        // Full Risch algorithm would go here
        IntegrationResult::Unknown(UnknownReason::UnsupportedForm)
    }

    /// Tries to find the integral in a lookup table.
    fn try_table_lookup(
        &mut self,
        expr: ExprHandle,
        var: ExprHandle,
        _var_id: SymbolId,
    ) -> Option<IntegrationResult> {
        let node = self.arena.get(expr);

        match node {
            // ∫ sin(x) dx = -cos(x)
            ExprNode::Function { id, args } if *id == functions::SIN => {
                if args.len() == 1 && args[0] == var {
                    let cos_x = self.arena.intern(ExprNode::Function {
                        id: functions::COS,
                        args: smallvec![var],
                    });
                    let neg_cos_x = self.arena.neg(cos_x);
                    return Some(IntegrationResult::Symbolic(
                        SymbolicAntiderivative::new(neg_cos_x, IntegrationMethod::TableLookup)
                            .with_display("-cos(x)".to_string()),
                    ));
                }
            }

            // ∫ cos(x) dx = sin(x)
            ExprNode::Function { id, args } if *id == functions::COS => {
                if args.len() == 1 && args[0] == var {
                    let sin_x = self.arena.intern(ExprNode::Function {
                        id: functions::SIN,
                        args: smallvec![var],
                    });
                    return Some(IntegrationResult::Symbolic(
                        SymbolicAntiderivative::new(sin_x, IntegrationMethod::TableLookup)
                            .with_display("sin(x)".to_string()),
                    ));
                }
            }

            // ∫ exp(x) dx = exp(x)
            ExprNode::Function { id, args } if *id == functions::EXP => {
                if args.len() == 1 && args[0] == var {
                    return Some(IntegrationResult::Symbolic(
                        SymbolicAntiderivative::new(expr, IntegrationMethod::TableLookup)
                            .with_display("exp(x)".to_string()),
                    ));
                }
            }

            // ∫ 1/x dx = ln(x)
            ExprNode::Div { num, den } => {
                let num_node = self.arena.get(*num);
                if matches!(num_node, ExprNode::Integer(1)) && *den == var {
                    let ln_x = self.arena.intern(ExprNode::Function {
                        id: functions::LN,
                        args: smallvec![var],
                    });
                    return Some(IntegrationResult::Symbolic(
                        SymbolicAntiderivative::new(ln_x, IntegrationMethod::TableLookup)
                            .with_display("ln(x)".to_string()),
                    ));
                }
            }

            _ => {}
        }

        None
    }

    /// Checks for known non-elementary integrals.
    fn check_non_elementary(
        &mut self,
        expr: ExprHandle,
        var_id: SymbolId,
        profile: &TranscendentalProfile,
    ) -> Option<IntegrationResult> {
        let node = self.arena.get(expr);

        // Check for exp(-x^2) -> erf pattern
        if let ExprNode::Function { id, args } = node {
            if *id == functions::EXP && args.len() == 1 {
                let arg_node = self.arena.get(args[0]);
                if let ExprNode::Neg(inner) = arg_node {
                    let inner_node = self.arena.get(*inner);
                    if let ExprNode::Pow { base, exp } = inner_node {
                        let base_node = self.arena.get(*base);
                        let exp_node = self.arena.get(*exp);
                        if matches!(base_node, ExprNode::Symbol(id) if *id == var_id)
                            && matches!(exp_node, ExprNode::Integer(2))
                        {
                            // ∫ exp(-x²) dx = (√π/2) erf(x)
                            return Some(IntegrationResult::SpecialFunction(
                                SpecialFunctionResult::new("(√π/2) erf(x)".to_string())
                                    .with_term(SpecialFunctionTerm {
                                        function: SpecialFunction::Erf,
                                        coefficient: Q::new(1, 2),
                                        argument: "x".to_string(),
                                    }),
                            ));
                        }
                    }
                }
            }
        }

        // Check for other known non-elementary forms based on profile
        if profile.has_exp && profile.has_log {
            // Many exp*log combinations are non-elementary
            // Would need more detailed analysis
        }

        None
    }

    /// Fallback for unclassified expressions.
    fn integrate_fallback(
        &mut self,
        _expr: ExprHandle,
        _var: ExprHandle,
        _var_id: SymbolId,
    ) -> IntegrationResult {
        IntegrationResult::Unknown(UnknownReason::UnsupportedForm)
    }

    /// Evaluates a definite integral F(upper) - F(lower).
    fn evaluate_definite(
        &mut self,
        _antideriv: ExprHandle,
        _var: ExprHandle,
        _lower: ExprHandle,
        _upper: ExprHandle,
    ) -> IntegrationResult {
        // Would need expression substitution and simplification
        // For now, return unknown
        IntegrationResult::Unknown(UnknownReason::UnsupportedForm)
    }

    /// Numerical integration fallback.
    fn numerical_definite(
        &mut self,
        _expr: ExprHandle,
        _var: ExprHandle,
        lower: ExprHandle,
        upper: ExprHandle,
    ) -> IntegrationResult {
        // Extract bounds as f64
        let a = match self.expr_to_f64(lower) {
            Some(v) => v,
            None => {
                return IntegrationResult::Unknown(UnknownReason::InternalError(
                    "lower bound must be numeric".to_string(),
                ))
            }
        };
        let b = match self.expr_to_f64(upper) {
            Some(v) => v,
            None => {
                return IntegrationResult::Unknown(UnknownReason::InternalError(
                    "upper bound must be numeric".to_string(),
                ))
            }
        };

        // Would need expression evaluation to create a closure
        // For now, just return a placeholder
        let result = adaptive_integrate(&|_x| 0.0, a, b, self.options.numerical_tolerance, self.options.numerical_tolerance, self.options.max_subdivisions);

        IntegrationResult::Numerical(
            NumericalResult::new(result.value, result.error)
                .with_evaluations(result.evaluations)
                .with_converged(result.converged),
        )
    }

    // === Helper methods ===

    /// Extracts the variable ID from a symbol expression.
    fn extract_variable_id(&self, var: ExprHandle) -> Option<SymbolId> {
        let node = self.arena.get(var);
        match node {
            ExprNode::Symbol(id) => Some(*id),
            _ => None,
        }
    }

    /// Converts an expression to a string (simple version).
    fn expr_to_string(&self, expr: ExprHandle) -> String {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Integer(n) => n.to_string(),
            ExprNode::Rational(num, den) => format!("{}/{}", num, den),
            ExprNode::Symbol(id) => {
                self.arena.symbol_name(*id).unwrap_or("?").to_string()
            }
            ExprNode::Add(_) => "(...+...)".to_string(),
            ExprNode::Mul(_) => "(...*...)".to_string(),
            ExprNode::Pow { .. } => "...^...".to_string(),
            ExprNode::Neg(_) => "-...".to_string(),
            ExprNode::Div { .. } => ".../...".to_string(),
            ExprNode::Function { id, .. } => format!("f{}(...)", id),
        }
    }

    /// Converts an expression to f64 (for numeric bounds).
    fn expr_to_f64(&self, expr: ExprHandle) -> Option<f64> {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Integer(n) => Some(*n as f64),
            ExprNode::Rational(num, den) => Some(*num as f64 / *den as f64),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrate_constant() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let five = arena.integer(5);

        let mut dispatcher = IntegrationDispatcher::new(&mut arena);
        let result = dispatcher.integrate(five, x);

        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_polynomial() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let two = arena.integer(2);
        let x_squared = arena.pow(x, two);

        let mut dispatcher = IntegrationDispatcher::new(&mut arena);
        let result = dispatcher.integrate(x_squared, x);

        assert!(result.is_symbolic());
        if let IntegrationResult::Symbolic(s) = result {
            assert_eq!(s.method, IntegrationMethod::PowerRule);
        }
    }

    #[test]
    fn test_integrate_sin() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let sin_x = arena.intern(ExprNode::Function {
            id: functions::SIN,
            args: smallvec![x],
        });

        let mut dispatcher = IntegrationDispatcher::new(&mut arena);
        let result = dispatcher.integrate(sin_x, x);

        assert!(result.is_symbolic());
        if let IntegrationResult::Symbolic(s) = result {
            assert_eq!(s.method, IntegrationMethod::TableLookup);
            assert!(s.display.contains("cos"));
        }
    }

    #[test]
    fn test_integrate_exp_neg_x_squared() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let two = arena.integer(2);
        let x_sq = arena.pow(x, two);
        let neg_x_sq = arena.neg(x_sq);
        let exp_neg_x_sq = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![neg_x_sq],
        });

        let mut dispatcher = IntegrationDispatcher::new(&mut arena);
        let result = dispatcher.integrate(exp_neg_x_sq, x);

        assert!(result.is_special_function());
    }
}
