//! The main simplification engine.
//!
//! This module provides the high-level API for simplifying expressions
//! using equality saturation.

use egg::{Extractor, RecExpr, Rewrite, Runner};

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
    /// Rewrite rules.
    rules: Vec<Rewrite<TertiusLang, ()>>,
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
            rules: rules::all_rules(),
        }
    }

    /// Creates a simplifier with custom configuration.
    #[must_use]
    pub fn with_config(config: SimplifierConfig) -> Self {
        Self {
            config,
            rules: rules::all_rules(),
        }
    }

    /// Sets custom rules (replaces default rules).
    #[must_use]
    pub fn with_rules(mut self, rules: Vec<Rewrite<TertiusLang, ()>>) -> Self {
        self.rules = rules;
        self
    }

    /// Adds rules to the existing set.
    pub fn add_rules(&mut self, rules: impl IntoIterator<Item = Rewrite<TertiusLang, ()>>) {
        self.rules.extend(rules);
    }

    /// Simplifies an expression given as a string.
    ///
    /// # Errors
    ///
    /// Returns an error if the expression cannot be parsed.
    pub fn simplify_str(&self, expr: &str) -> Result<String, String> {
        let parsed: RecExpr<TertiusLang> = expr
            .parse()
            .map_err(|e| format!("parse error: {e}"))?;

        let simplified = self.simplify(&parsed);
        Ok(simplified.to_string())
    }

    /// Simplifies a parsed expression.
    #[must_use]
    pub fn simplify(&self, expr: &RecExpr<TertiusLang>) -> RecExpr<TertiusLang> {
        let runner = Runner::default()
            .with_expr(expr)
            .with_iter_limit(self.config.iter_limit)
            .with_node_limit(self.config.node_limit)
            .with_time_limit(std::time::Duration::from_secs(self.config.time_limit_secs))
            .run(&self.rules);

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
        let runner = Runner::default()
            .with_expr(expr)
            .with_iter_limit(self.config.iter_limit)
            .with_node_limit(self.config.node_limit)
            .with_time_limit(std::time::Duration::from_secs(self.config.time_limit_secs))
            .run(&self.rules);

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

#[cfg(test)]
mod tests {
    use super::*;

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
        let result = simplifier.simplify_str("(+ (^ (sin x) 2) (^ (cos x) 2))").unwrap();
        assert_eq!(result, "1");
    }

    #[test]
    fn test_simplify_exp_log() {
        let simplifier = Simplifier::new();

        // exp(ln(x)) = x
        assert_eq!(simplifier.simplify_str("(exp (ln x))").unwrap(), "x");

        // ln(exp(x)) = x
        assert_eq!(simplifier.simplify_str("(ln (exp x))").unwrap(), "x");
    }
}
