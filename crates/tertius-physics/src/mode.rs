//! High-level physics-mode solving API.

use tertius_core::assumptions::{Assumption, AssumptionSet};
use tertius_solve::solve_first_order_linear_constant;
use tertius_units::Unit;

use crate::trace::DerivationTrace;

/// Inputs for solving `y' + a y = b`.
#[derive(Clone, Debug)]
pub struct FirstOrderLinearRequest {
    /// Constant `a`.
    pub a: f64,
    /// Constant `b`.
    pub b: f64,
    /// Optional initial value `y(0)`.
    pub y0: Option<f64>,
    /// Unit for independent variable `x`.
    pub x_unit: Unit,
    /// Unit for dependent variable `y`.
    pub y_unit: Unit,
    /// Unit of coefficient `a`.
    pub a_unit: Unit,
    /// Unit of coefficient `b`.
    pub b_unit: Unit,
}

/// Physics-mode solve result with UX-centric payload.
#[derive(Clone, Debug)]
pub struct PhysicsSolveResult {
    /// Human-readable final expression.
    pub expression: String,
    /// LaTeX representation of the final expression.
    pub expression_latex: String,
    /// Diagnostics (unit mismatch, unsafe assumptions, etc.).
    pub diagnostics: Vec<String>,
    /// Rule-level derivation trace.
    pub derivation: DerivationTrace,
}

/// Physics-mode facade that combines assumptions, units, and derivation traces.
#[derive(Clone, Debug, Default)]
pub struct PhysicsMode {
    assumptions: AssumptionSet,
}

impl PhysicsMode {
    /// Creates a default physics mode.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an assumption.
    pub fn assume(&mut self, symbol: &str, assumption: Assumption) {
        self.assumptions.assume(symbol, assumption);
    }

    /// Solves `y' + a y = b` and returns a physics-facing result bundle.
    #[must_use]
    pub fn solve_first_order_linear(
        &self,
        request: &FirstOrderLinearRequest,
    ) -> PhysicsSolveResult {
        let solution = solve_first_order_linear_constant(request.a, request.b, request.y0);
        let _sample = solution.eval(0.0);

        let mut diagnostics = Vec::new();
        let expected_a = Unit::new(
            "1/x",
            1.0 / request.x_unit.scale_to_si,
            request.x_unit.dimension.powi(-1),
        );
        if !request.a_unit.is_compatible_with(&expected_a) {
            diagnostics.push(format!(
                "unit mismatch: coefficient a has {}, expected inverse {}",
                request.a_unit.symbol, request.x_unit.symbol
            ));
        }
        let expected_b = request.y_unit.div(&request.x_unit, "y/x");
        if !request.b_unit.is_compatible_with(&expected_b) {
            diagnostics.push(format!(
                "unit mismatch: coefficient b has {}, expected {}/{}",
                request.b_unit.symbol, request.y_unit.symbol, request.x_unit.symbol
            ));
        }

        let (expression, expression_latex) =
            build_first_order_solution_strings(request.a, request.b, request.y0);

        let mut trace = DerivationTrace::new("First-Order Linear ODE");
        trace.push_step(
            "Normalize equation",
            "y' + a y = b",
            "y' = b - a y",
            Vec::new(),
        );
        let mut cond = Vec::new();
        if self.assumptions.has("x", Assumption::Real)
            || self.assumptions.has("x", Assumption::Positive)
        {
            cond.push("x is real".to_string());
        } else {
            cond.push("assume x is real".to_string());
        }
        trace.push_step("Integrating factor", "mu = e^{∫a dx}", "mu = e^{a x}", cond);
        trace.push_step(
            "Exact derivative form",
            "y' + a y = b",
            "(mu y)' = b mu",
            Vec::new(),
        );
        if let Some(y0) = request.y0 {
            trace.push_step(
                "Apply initial condition",
                &format!("y(0) = {}", fmt_num(y0)),
                "solve for integration constant",
                Vec::new(),
            );
        }
        trace.final_result = expression.clone();

        PhysicsSolveResult {
            expression,
            expression_latex,
            diagnostics,
            derivation: trace,
        }
    }
}

fn build_first_order_solution_strings(a: f64, b: f64, y0: Option<f64>) -> (String, String) {
    if a.abs() < 1e-12 {
        return match y0 {
            Some(y_init) => (
                format!("y(x) = {} x + {}", fmt_num(b), fmt_num(y_init)),
                format!("y(x) = {} x + {}", fmt_num(b), fmt_num(y_init)),
            ),
            None => (
                format!("y(x) = {} x + C", fmt_num(b)),
                format!("y(x) = {} x + C", fmt_num(b)),
            ),
        };
    }

    let steady = b / a;
    match y0 {
        Some(y_init) => {
            let c = y_init - steady;
            (
                format!(
                    "y(x) = {} + {} e^(-{} x)",
                    fmt_num(steady),
                    fmt_num(c),
                    fmt_num(a)
                ),
                format!(
                    "y(x) = {} + {} e^{{-{} x}}",
                    fmt_num(steady),
                    fmt_num(c),
                    fmt_num(a)
                ),
            )
        }
        None => (
            format!("y(x) = {} + C e^(-{} x)", fmt_num(steady), fmt_num(a)),
            format!("y(x) = {} + C e^{{-{} x}}", fmt_num(steady), fmt_num(a)),
        ),
    }
}

fn fmt_num(v: f64) -> String {
    let rounded = v.round();
    if (v - rounded).abs() < 1e-12 {
        format!("{}", rounded as i64)
    } else {
        let mut s = format!("{v:.6}");
        while s.contains('.') && s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_request() -> FirstOrderLinearRequest {
        FirstOrderLinearRequest {
            a: 2.0,
            b: 4.0,
            y0: Some(1.0),
            x_unit: Unit::second(),
            y_unit: Unit::meter(),
            a_unit: Unit::new("1/s", 1.0, Unit::second().dimension.powi(-1)),
            b_unit: Unit::new(
                "m/s",
                1.0,
                Unit::meter().div(&Unit::second(), "m/s").dimension,
            ),
        }
    }

    #[test]
    fn test_solve_first_order_linear_produces_derivation() {
        let mut mode = PhysicsMode::new();
        mode.assume("x", Assumption::Real);
        let result = mode.solve_first_order_linear(&base_request());
        assert!(result.expression.contains("y(x)"));
        assert!(result.expression_latex.contains("y(x)"));
        assert!(!result.derivation.steps.is_empty());
        assert!(result
            .derivation
            .steps
            .iter()
            .any(|s| s.rule_applied.contains("Integrating factor")));
    }

    #[test]
    fn test_unit_mismatch_emits_diagnostic() {
        let mode = PhysicsMode::new();
        let mut req = base_request();
        req.a_unit = Unit::meter();
        let result = mode.solve_first_order_linear(&req);
        assert!(!result.diagnostics.is_empty());
        assert!(result.diagnostics.iter().any(|d| d.contains("unit")));
    }
}
