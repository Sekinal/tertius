//! Most Rapidly Varying (MRV) subexpression detection.
//!
//! The MRV set contains subexpressions that grow fastest as x → ∞.
//! This is the key insight of the Gruntz algorithm: by rewriting in terms
//! of the fastest-growing subexpression, we can reduce the problem.

use crate::comparison::{ComparisonClass, GrowthRate};
use std::collections::HashSet;
use tertius_core::{ExprArena, ExprHandle, ExprNode};
use tertius_core::expr::functions;

/// A set of Most Rapidly Varying subexpressions.
#[derive(Clone, Debug)]
pub struct MrvSet {
    /// The expressions in the MRV set.
    exprs: Vec<ExprHandle>,
    /// The common growth rate of all expressions in the set.
    growth_rate: GrowthRate,
}

impl MrvSet {
    /// Creates an empty MRV set.
    pub fn empty() -> Self {
        Self {
            exprs: Vec::new(),
            growth_rate: GrowthRate::Constant,
        }
    }

    /// Creates a singleton MRV set.
    pub fn singleton(expr: ExprHandle, rate: GrowthRate) -> Self {
        Self {
            exprs: vec![expr],
            growth_rate: rate,
        }
    }

    /// Returns true if the set is empty.
    pub fn is_empty(&self) -> bool {
        self.exprs.is_empty()
    }

    /// Returns the expressions in the MRV set.
    pub fn exprs(&self) -> &[ExprHandle] {
        &self.exprs
    }

    /// Returns the growth rate.
    pub fn growth_rate(&self) -> &GrowthRate {
        &self.growth_rate
    }

    /// Merges two MRV sets based on their growth rates.
    ///
    /// Returns the set of faster-growing expressions.
    /// If both have the same growth rate, returns their union.
    pub fn merge(self, other: Self) -> Self {
        if self.is_empty() {
            return other;
        }
        if other.is_empty() {
            return self;
        }

        match self.growth_rate.compare(&other.growth_rate) {
            ComparisonClass::LessThan => other,
            ComparisonClass::GreaterThan => self,
            ComparisonClass::Comparable => {
                // Same growth rate: merge the sets
                let mut combined = self.exprs;
                let seen: HashSet<ExprHandle> = combined.iter().copied().collect();
                for e in other.exprs {
                    if !seen.contains(&e) {
                        combined.push(e);
                    }
                }
                Self {
                    exprs: combined,
                    growth_rate: self.growth_rate,
                }
            }
            ComparisonClass::Unknown => {
                // Unknown: keep both (conservative)
                let mut combined = self.exprs;
                combined.extend(other.exprs);
                Self {
                    exprs: combined,
                    growth_rate: GrowthRate::Unknown,
                }
            }
        }
    }

    /// Returns the dominant expression (first in the set).
    pub fn dominant(&self) -> Option<ExprHandle> {
        self.exprs.first().copied()
    }
}

/// Context for MRV computation.
pub struct MrvContext<'a> {
    arena: &'a ExprArena,
    variable: ExprHandle,
}

impl<'a> MrvContext<'a> {
    /// Creates a new MRV context.
    pub fn new(arena: &'a ExprArena, variable: ExprHandle) -> Self {
        Self { arena, variable }
    }

    /// Determines the growth rate of an expression.
    pub fn growth_rate(&self, expr: ExprHandle) -> GrowthRate {
        let e = self.arena.get(expr);
        match e {
            ExprNode::Integer(_) | ExprNode::Rational(_, _) => GrowthRate::Constant,

            ExprNode::Symbol(_) => {
                // Check if this is our variable
                if expr == self.variable {
                    GrowthRate::Polynomial(1)
                } else {
                    // Other symbols are treated as constants
                    GrowthRate::Constant
                }
            }

            ExprNode::Add(args) | ExprNode::Mul(args) => {
                // For sum/product, the growth rate is determined by the fastest term
                let mut max_rate = GrowthRate::Constant;
                for &arg in args {
                    let rate = self.growth_rate(arg);
                    if rate.compare(&max_rate) == ComparisonClass::GreaterThan {
                        max_rate = rate;
                    }
                }
                max_rate
            }

            ExprNode::Pow { base, exp } => {
                let base_rate = self.growth_rate(*base);
                let exp_rate = self.growth_rate(*exp);

                // If exponent has exponential growth, result is super-exponential
                if matches!(
                    exp_rate,
                    GrowthRate::Exponential { .. } | GrowthRate::SuperExponential
                ) {
                    return GrowthRate::SuperExponential;
                }

                // If base is x and exponent is constant, it's polynomial
                if *base == self.variable {
                    if let ExprNode::Integer(n) = self.arena.get(*exp) {
                        if let Ok(n_i32) = i32::try_from(*n) {
                            return GrowthRate::Polynomial(n_i32);
                        }
                    }
                }

                // If exponent contains x, it's exponential
                if self.contains_variable(*exp) {
                    GrowthRate::Exponential {
                        coeff: 1.0,
                        power: 1,
                    }
                } else {
                    base_rate
                }
            }

            ExprNode::Function { id, args } => {
                match *id {
                    functions::EXP => {
                        if let Some(&arg) = args.first() {
                            let arg_rate = self.growth_rate(arg);
                            match arg_rate {
                                GrowthRate::Constant => GrowthRate::Constant,
                                GrowthRate::Logarithmic(_) => GrowthRate::Polynomial(1),
                                GrowthRate::Polynomial(n) => GrowthRate::Exponential {
                                    coeff: 1.0,
                                    power: n,
                                },
                                GrowthRate::Exponential { .. } | GrowthRate::SuperExponential => {
                                    GrowthRate::SuperExponential
                                }
                                GrowthRate::Unknown => GrowthRate::Unknown,
                            }
                        } else {
                            GrowthRate::Unknown
                        }
                    }

                    functions::LN => {
                        if let Some(&arg) = args.first() {
                            let arg_rate = self.growth_rate(arg);
                            match arg_rate {
                                GrowthRate::Constant => GrowthRate::Constant,
                                GrowthRate::Polynomial(n) if n > 0 => GrowthRate::Logarithmic(1),
                                GrowthRate::Exponential { .. } => GrowthRate::Polynomial(1),
                                GrowthRate::SuperExponential => GrowthRate::Exponential {
                                    coeff: 1.0,
                                    power: 1,
                                },
                                _ => GrowthRate::Unknown,
                            }
                        } else {
                            GrowthRate::Unknown
                        }
                    }

                    functions::SIN | functions::COS => {
                        // Bounded functions
                        GrowthRate::Constant
                    }

                    _ => GrowthRate::Unknown,
                }
            }

            ExprNode::Neg(arg) => self.growth_rate(*arg),

            ExprNode::Div { num, den } => {
                let num_rate = self.growth_rate(*num);
                let den_rate = self.growth_rate(*den);
                // For division, compare rates
                match num_rate.compare(&den_rate) {
                    ComparisonClass::GreaterThan => num_rate,
                    ComparisonClass::LessThan => {
                        // Decays to zero, but rate depends on denominator
                        GrowthRate::Constant
                    }
                    ComparisonClass::Comparable => GrowthRate::Constant, // Same order → constant
                    ComparisonClass::Unknown => GrowthRate::Unknown,
                }
            }
        }
    }

    /// Checks if an expression contains the limit variable.
    pub fn contains_variable(&self, expr: ExprHandle) -> bool {
        if expr == self.variable {
            return true;
        }

        let e = self.arena.get(expr);
        match e {
            ExprNode::Integer(_) | ExprNode::Rational(_, _) | ExprNode::Symbol(_) => false,
            ExprNode::Add(args) | ExprNode::Mul(args) => {
                args.iter().any(|&a| self.contains_variable(a))
            }
            ExprNode::Pow { base, exp } => {
                self.contains_variable(*base) || self.contains_variable(*exp)
            }
            ExprNode::Function { args, .. } => args.iter().any(|&a| self.contains_variable(a)),
            ExprNode::Neg(arg) => self.contains_variable(*arg),
            ExprNode::Div { num, den } => {
                self.contains_variable(*num) || self.contains_variable(*den)
            }
        }
    }
}

/// Finds the MRV set of an expression as x → ∞.
///
/// The MRV set contains the subexpressions that grow fastest.
/// All expressions in the set have the same asymptotic growth rate.
pub fn find_mrv(arena: &ExprArena, expr: ExprHandle, variable: ExprHandle) -> MrvSet {
    let ctx = MrvContext::new(arena, variable);
    find_mrv_recursive(&ctx, expr)
}

/// Recursive MRV computation.
fn find_mrv_recursive(ctx: &MrvContext<'_>, expr: ExprHandle) -> MrvSet {
    let e = ctx.arena.get(expr);

    match e {
        ExprNode::Integer(_) | ExprNode::Rational(_, _) => MrvSet::empty(),

        ExprNode::Symbol(_) => {
            if expr == ctx.variable {
                // x has polynomial(1) growth
                MrvSet::singleton(expr, GrowthRate::Polynomial(1))
            } else {
                MrvSet::empty()
            }
        }

        ExprNode::Add(args) | ExprNode::Mul(args) => {
            let mut result = MrvSet::empty();
            for &arg in args {
                let arg_mrv = find_mrv_recursive(ctx, arg);
                result = result.merge(arg_mrv);
            }
            result
        }

        ExprNode::Pow { base, exp } => {
            let base_mrv = find_mrv_recursive(ctx, *base);
            let exp_mrv = find_mrv_recursive(ctx, *exp);

            // If exp contains x, the whole expression might be in MRV
            if ctx.contains_variable(*exp) {
                let rate = ctx.growth_rate(expr);
                let expr_mrv = MrvSet::singleton(expr, rate);
                base_mrv.merge(exp_mrv).merge(expr_mrv)
            } else {
                base_mrv.merge(exp_mrv)
            }
        }

        ExprNode::Function { id, args } => {
            // Get MRV from arguments first
            let mut result = MrvSet::empty();
            for &arg in args {
                let arg_mrv = find_mrv_recursive(ctx, arg);
                result = result.merge(arg_mrv);
            }

            match *id {
                functions::EXP => {
                    // exp(f(x)) is in MRV if f(x) → ±∞
                    if let Some(&arg) = args.first() {
                        let arg_rate = ctx.growth_rate(arg);
                        if matches!(
                            arg_rate,
                            GrowthRate::Polynomial(_) | GrowthRate::Exponential { .. }
                        ) {
                            let rate = ctx.growth_rate(expr);
                            let expr_mrv = MrvSet::singleton(expr, rate);
                            result = result.merge(expr_mrv);
                        }
                    }
                }
                _ => {}
            }

            result
        }

        ExprNode::Neg(arg) => find_mrv_recursive(ctx, *arg),

        ExprNode::Div { num, den } => {
            let num_mrv = find_mrv_recursive(ctx, *num);
            let den_mrv = find_mrv_recursive(ctx, *den);
            num_mrv.merge(den_mrv)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;

    fn setup_arena() -> (ExprArena, ExprHandle) {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        (arena, x)
    }

    #[test]
    fn test_mrv_constant() {
        let (mut arena, x) = setup_arena();
        let c = arena.integer(5);

        let mrv = find_mrv(&arena, c, x);
        assert!(mrv.is_empty());
    }

    #[test]
    fn test_mrv_variable() {
        let (arena, x) = setup_arena();

        let mrv = find_mrv(&arena, x, x);
        assert_eq!(mrv.exprs().len(), 1);
        assert_eq!(mrv.exprs()[0], x);
    }

    #[test]
    fn test_mrv_polynomial() {
        let (mut arena, x) = setup_arena();
        // x^2
        let two = arena.integer(2);
        let x_squared = arena.pow(x, two);

        let mrv = find_mrv(&arena, x_squared, x);
        // MRV should contain x (the fastest growing subexpr)
        assert!(!mrv.is_empty());
    }

    #[test]
    fn test_mrv_exp() {
        let (mut arena, x) = setup_arena();
        // exp(x)
        let exp_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![x],
        });

        let mrv = find_mrv(&arena, exp_x, x);
        // MRV should contain exp(x)
        assert!(!mrv.is_empty());
        assert!(mrv.exprs().contains(&exp_x));
    }

    #[test]
    fn test_mrv_merge() {
        let mrv1 = MrvSet::singleton(ExprHandle::new(1), GrowthRate::Polynomial(2));
        let mrv2 = MrvSet::singleton(ExprHandle::new(2), GrowthRate::Exponential {
            coeff: 1.0,
            power: 1,
        });

        let merged = mrv1.merge(mrv2);
        // Exponential > Polynomial, so exp should dominate
        assert!(matches!(
            merged.growth_rate(),
            GrowthRate::Exponential { .. }
        ));
    }

    #[test]
    fn test_growth_rate_log() {
        let (mut arena, x) = setup_arena();
        let log_x = arena.intern(ExprNode::Function {
            id: functions::LN,
            args: smallvec![x],
        });

        let ctx = MrvContext::new(&arena, x);
        let rate = ctx.growth_rate(log_x);

        assert!(matches!(rate, GrowthRate::Logarithmic(1)));
    }

    #[test]
    fn test_growth_rate_hierarchy() {
        // log(x) < x < x^2 < exp(x) < exp(exp(x))
        let (mut arena, x) = setup_arena();

        let log_x = arena.intern(ExprNode::Function {
            id: functions::LN,
            args: smallvec![x],
        });
        let two = arena.integer(2);
        let x_squared = arena.pow(x, two);
        let exp_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![x],
        });
        let exp_exp_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![exp_x],
        });

        let ctx = MrvContext::new(&arena, x);

        let r_log = ctx.growth_rate(log_x);
        let r_x = ctx.growth_rate(x);
        let r_x2 = ctx.growth_rate(x_squared);
        let r_exp = ctx.growth_rate(exp_x);
        let r_exp_exp = ctx.growth_rate(exp_exp_x);

        assert_eq!(r_log.compare(&r_x), ComparisonClass::LessThan);
        assert_eq!(r_x.compare(&r_x2), ComparisonClass::LessThan);
        assert_eq!(r_x2.compare(&r_exp), ComparisonClass::LessThan);
        assert_eq!(r_exp.compare(&r_exp_exp), ComparisonClass::LessThan);
    }
}
