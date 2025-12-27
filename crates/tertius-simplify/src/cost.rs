//! Cost functions for expression extraction.
//!
//! After equality saturation, we need to pick the "best" expression
//! from each equivalence class. Cost functions define what "best" means.

use egg::{CostFunction, Id, Language};

use crate::language::TertiusLang;

/// A cost function that minimizes AST size.
///
/// This is the simplest cost function: prefer smaller expressions.
#[derive(Default)]
pub struct AstSizeCost;

impl CostFunction<TertiusLang> for AstSizeCost {
    type Cost = usize;

    fn cost<C>(&mut self, enode: &TertiusLang, mut costs: C) -> Self::Cost
    where
        C: FnMut(Id) -> Self::Cost,
    {
        let base_cost = match enode {
            // Atoms have cost 1
            TertiusLang::Num(_) | TertiusLang::Symbol(_) => 1,
            // Other nodes have cost 1 + children
            _ => 1,
        };

        enode.fold(base_cost, |sum, id| sum + costs(id))
    }
}

/// A cost function that prefers certain forms.
///
/// This can be customized to prefer:
/// - Factored forms
/// - Expanded forms
/// - Certain function representations
#[derive(Default)]
pub struct WeightedCost {
    /// Penalty for expanded products.
    pub expansion_penalty: usize,
    /// Bonus for factored forms.
    pub factoring_bonus: usize,
}

impl CostFunction<TertiusLang> for WeightedCost {
    type Cost = usize;

    fn cost<C>(&mut self, enode: &TertiusLang, mut costs: C) -> Self::Cost
    where
        C: FnMut(Id) -> Self::Cost,
    {
        let base_cost = match enode {
            // Numbers are free
            TertiusLang::Num(_) => 0,
            // Symbols have cost 1
            TertiusLang::Symbol(_) => 1,
            // Basic ops
            TertiusLang::Add(_) | TertiusLang::Sub(_) => 2,
            TertiusLang::Mul(_) => 2,
            TertiusLang::Div(_) => 3,
            TertiusLang::Neg(_) => 1,
            TertiusLang::Pow(_) => 3,
            // Transcendental functions are expensive
            TertiusLang::Sin(_)
            | TertiusLang::Cos(_)
            | TertiusLang::Tan(_)
            | TertiusLang::Exp(_)
            | TertiusLang::Ln(_) => 5,
            // Others
            _ => 3,
        };

        enode.fold(base_cost, |sum, id| sum + costs(id))
    }
}

/// A cost function optimized for numerical evaluation.
///
/// Prefers forms that minimize floating-point operations.
#[derive(Default)]
pub struct NumericalCost;

impl CostFunction<TertiusLang> for NumericalCost {
    type Cost = usize;

    fn cost<C>(&mut self, enode: &TertiusLang, mut costs: C) -> Self::Cost
    where
        C: FnMut(Id) -> Self::Cost,
    {
        let base_cost = match enode {
            TertiusLang::Num(_) => 0,
            TertiusLang::Symbol(_) => 0,
            // Prefer fused multiply-add forms
            TertiusLang::Add(_) => 1,
            TertiusLang::Mul(_) => 1,
            TertiusLang::Sub(_) => 1,
            TertiusLang::Neg(_) => 0, // Free on most architectures
            // Division is expensive
            TertiusLang::Div(_) => 10,
            // Powers: depends on exponent
            TertiusLang::Pow(_) => 5,
            // Transcendental: very expensive
            TertiusLang::Sin(_) | TertiusLang::Cos(_) | TertiusLang::Tan(_) => 20,
            TertiusLang::Exp(_) | TertiusLang::Ln(_) => 15,
            TertiusLang::Sqrt(_) => 5,
            _ => 5,
        };

        enode.fold(base_cost, |sum, id| sum + costs(id))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egg::{Extractor, RecExpr, Runner};

    #[test]
    fn test_ast_size() {
        let expr: RecExpr<TertiusLang> = "(+ x 0)".parse().unwrap();
        let rules = crate::rules::arithmetic::rules();

        let runner = Runner::default().with_expr(&expr).run(&rules);
        let extractor = Extractor::new(&runner.egraph, AstSizeCost);
        let (cost, best) = extractor.find_best(runner.roots[0]);

        // "x" should be cheaper than "(+ x 0)"
        assert_eq!(best.to_string(), "x");
        assert_eq!(cost, 1);
    }
}
