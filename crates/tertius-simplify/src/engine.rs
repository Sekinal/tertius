//! The main simplification engine.
//!
//! This module provides the high-level API for simplifying expressions
//! using equality saturation.

use egg::{Extractor, RecExpr, Rewrite, Runner};
use tertius_core::assumptions::AssumptionSet;

use crate::cost::AstSizeCost;
use crate::language::TertiusLang;
use crate::rules;

/// Configuration for the simplification engine.
#[derive(Clone, Debug)]
pub struct SimplifierConfig {
    /// Maximum number of iterations.
    pub iter_limit: usize,
    /// Maximum number of nodes in the e-graph.
    pub node_limit: usize,
    /// Time limit in seconds.
    pub time_limit_secs: u64,
}

impl Default for SimplifierConfig {
    fn default() -> Self {
        Self {
            iter_limit: 30,
            node_limit: 100_000,
            time_limit_secs: 10,
        }
    }
}

/// The main simplification engine.
#[derive(Clone)]
pub struct Simplifier {
    /// Configuration.
    config: SimplifierConfig,
    /// Rewrite rules that are always branch-safe.
    base_rules: Vec<Rewrite<TertiusLang, ()>>,
    /// Rewrite rules requiring explicit assumptions.
    real_assumption_rules: Vec<Rewrite<TertiusLang, ()>>,
}

impl Default for Simplifier {
    fn default() -> Self {
        Self::new()
    }
}

impl Simplifier {
    /// Creates a new simplifier with default rules.
    #[must_use]
    pub fn new() -> Self {
        Self {
            config: SimplifierConfig::default(),
            base_rules: rules::all_rules(),
            real_assumption_rules: rules::real_assumption_rules(),
        }
    }

    /// Creates a simplifier with custom configuration.
    #[must_use]
    pub fn with_config(config: SimplifierConfig) -> Self {
        Self {
            config,
            base_rules: rules::all_rules(),
            real_assumption_rules: rules::real_assumption_rules(),
        }
    }

    /// Sets custom rules (replaces default rules).
    #[must_use]
    pub fn with_rules(mut self, rules: Vec<Rewrite<TertiusLang, ()>>) -> Self {
        self.base_rules = rules;
        self
    }

    /// Adds rules to the existing set.
    pub fn add_rules(&mut self, rules: impl IntoIterator<Item = Rewrite<TertiusLang, ()>>) {
        self.base_rules.extend(rules);
    }

    /// Simplifies an expression given as a string.
    ///
    /// # Errors
    ///
    /// Returns an error if the expression cannot be parsed.
    pub fn simplify_str(&self, expr: &str) -> Result<String, String> {
        let parsed: RecExpr<TertiusLang> = expr.parse().map_err(|e| format!("parse error: {e}"))?;

        let simplified = self.simplify(&parsed);
        Ok(simplified.to_string())
    }

    /// Simplifies a string expression with explicit symbol assumptions.
    ///
    /// Domain-sensitive rewrites are enabled only when all symbols in the
    /// expression are assumed real (or positive).
    pub fn simplify_str_with_assumptions(
        &self,
        expr: &str,
        assumptions: &AssumptionSet,
    ) -> Result<String, String> {
        let parsed: RecExpr<TertiusLang> = expr.parse().map_err(|e| format!("parse error: {e}"))?;

        let simplified = self.simplify_with_assumptions(&parsed, Some(assumptions));
        Ok(simplified.to_string())
    }

    /// Simplifies a parsed expression.
    #[must_use]
    pub fn simplify(&self, expr: &RecExpr<TertiusLang>) -> RecExpr<TertiusLang> {
        self.simplify_with_assumptions(expr, None)
    }

    /// Simplifies with optional assumptions for domain-sensitive rewrites.
    #[must_use]
    pub fn simplify_with_assumptions(
        &self,
        expr: &RecExpr<TertiusLang>,
        assumptions: Option<&AssumptionSet>,
    ) -> RecExpr<TertiusLang> {
        let rules = self.rules_for_expr(expr, assumptions);
        let runner = Runner::default()
            .with_expr(expr)
            .with_iter_limit(self.config.iter_limit)
            .with_node_limit(self.config.node_limit)
            .with_time_limit(std::time::Duration::from_secs(self.config.time_limit_secs))
            .run(&rules);

        let extractor = Extractor::new(&runner.egraph, AstSizeCost);
        let (_, best) = extractor.find_best(runner.roots[0]);
        best
    }

    /// Simplifies and returns both the result and statistics.
    #[must_use]
    pub fn simplify_with_stats(
        &self,
        expr: &RecExpr<TertiusLang>,
    ) -> (RecExpr<TertiusLang>, SimplificationStats) {
        let rules = self.rules_for_expr(expr, None);
        let runner = Runner::default()
            .with_expr(expr)
            .with_iter_limit(self.config.iter_limit)
            .with_node_limit(self.config.node_limit)
            .with_time_limit(std::time::Duration::from_secs(self.config.time_limit_secs))
            .run(&rules);

        let stats = SimplificationStats {
            iterations: runner.iterations.len(),
            egraph_nodes: runner.egraph.total_number_of_nodes(),
            egraph_classes: runner.egraph.number_of_classes(),
            stop_reason: format!("{:?}", runner.stop_reason),
        };

        let extractor = Extractor::new(&runner.egraph, AstSizeCost);
        let (_, best) = extractor.find_best(runner.roots[0]);

        (best, stats)
    }

    fn rules_for_expr(
        &self,
        expr: &RecExpr<TertiusLang>,
        assumptions: Option<&AssumptionSet>,
    ) -> Vec<Rewrite<TertiusLang, ()>> {
        let mut rules = self.base_rules.clone();
        if let Some(assumptions) = assumptions {
            let symbol_names = collect_symbols(expr);
            if assumptions.all_real(symbol_names.iter().map(std::string::String::as_str)) {
                rules.extend(self.real_assumption_rules.iter().cloned());
            }
        }
        rules
    }
}

/// Statistics about the simplification process.
#[derive(Clone, Debug)]
pub struct SimplificationStats {
    /// Number of iterations run.
    pub iterations: usize,
    /// Total nodes in the e-graph.
    pub egraph_nodes: usize,
    /// Number of equivalence classes.
    pub egraph_classes: usize,
    /// Reason the runner stopped.
    pub stop_reason: String,
}

fn collect_symbols(expr: &RecExpr<TertiusLang>) -> Vec<String> {
    let mut out = Vec::new();
    for node in expr.as_ref() {
        if let TertiusLang::Symbol(sym) = node {
            out.push(sym.to_string());
        }
    }
    out.sort();
    out.dedup();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_core::assumptions::{Assumption, AssumptionSet};

    #[test]
    fn test_simplify_basic() {
        let simplifier = Simplifier::new();

        // x + 0 = x
        assert_eq!(simplifier.simplify_str("(+ x 0)").unwrap(), "x");

        // x * 1 = x
        assert_eq!(simplifier.simplify_str("(* x 1)").unwrap(), "x");

        // x * 0 = 0
        assert_eq!(simplifier.simplify_str("(* x 0)").unwrap(), "0");
    }

    #[test]
    fn test_simplify_trig() {
        let simplifier = Simplifier::new();

        // sin²(x) + cos²(x) = 1
        let result = simplifier
            .simplify_str("(+ (^ (sin x) 2) (^ (cos x) 2))")
            .unwrap();
        assert_eq!(result, "1");
    }

    #[test]
    fn test_simplify_exp_log() {
        let simplifier = Simplifier::new();

        // exp(ln(x)) = x
        assert_eq!(simplifier.simplify_str("(exp (ln x))").unwrap(), "x");

        // ln(exp(x)) is branch-sensitive and should not simplify by default
        assert_eq!(
            simplifier.simplify_str("(ln (exp x))").unwrap(),
            "(ln (exp x))"
        );
    }

    #[test]
    fn test_simplify_ln_exp_with_real_assumption() {
        let simplifier = Simplifier::new();
        let mut assumptions = AssumptionSet::new();
        assumptions.assume("x", Assumption::Real);

        let result = simplifier
            .simplify_str_with_assumptions("(ln (exp x))", &assumptions)
            .unwrap();
        assert_eq!(result, "x");
    }
}
