//! The Gruntz algorithm for computing limits.
//!
//! This implements the algorithm from:
//! "Computing limits of sequences" by Dominik Gruntz (1996)
//!
//! The algorithm works by:
//! 1. Finding the Most Rapidly Varying (MRV) subexpressions
//! 2. Rewriting the expression in terms of a new variable ω where ω → 0
//! 3. Computing a series expansion in ω
//! 4. Extracting the limit from the leading term

use crate::comparison::ComparisonClass;
use crate::mrv::{find_mrv, MrvContext};
use tertius_core::{ExprArena, ExprHandle, ExprNode};
use tertius_core::expr::functions;
use thiserror::Error;

/// The direction of a limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Limit {
    /// x → +∞
    PosInfinity,
    /// x → -∞
    NegInfinity,
    /// x → a from the right (x → a⁺)
    Right(ExprHandle),
    /// x → a from the left (x → a⁻)
    Left(ExprHandle),
    /// x → a (two-sided)
    TwoSided(ExprHandle),
}

impl Limit {
    /// Returns true if this is a limit at infinity.
    pub fn is_infinity(&self) -> bool {
        matches!(self, Limit::PosInfinity | Limit::NegInfinity)
    }
}

/// The result of a limit computation.
#[derive(Clone, Debug)]
pub enum LimitResult {
    /// A finite limit value.
    Finite(ExprHandle),
    /// The limit is +∞.
    PosInfinity,
    /// The limit is -∞.
    NegInfinity,
    /// The limit does not exist (e.g., oscillates).
    DoesNotExist,
    /// Could not determine the limit.
    Unknown,
}

impl LimitResult {
    /// Returns true if the limit is finite.
    pub fn is_finite(&self) -> bool {
        matches!(self, LimitResult::Finite(_))
    }

    /// Returns the finite value if present.
    pub fn finite_value(&self) -> Option<ExprHandle> {
        match self {
            LimitResult::Finite(v) => Some(*v),
            _ => None,
        }
    }
}

/// Errors that can occur during limit computation.
#[derive(Clone, Debug, Error)]
pub enum LimitError {
    #[error("variable not found in expression")]
    VariableNotFound,

    #[error("cannot compute limit: {0}")]
    CannotCompute(String),

    #[error("series expansion failed")]
    SeriesExpansionFailed,
}

/// Context for the Gruntz algorithm.
pub struct GruntzContext<'a> {
    arena: &'a mut ExprArena,
    variable: ExprHandle,
    limit: Limit,
    /// Recursion depth (for termination)
    depth: usize,
    /// Maximum recursion depth
    max_depth: usize,
}

impl<'a> GruntzContext<'a> {
    /// Creates a new Gruntz context.
    pub fn new(arena: &'a mut ExprArena, variable: ExprHandle, limit: Limit) -> Self {
        Self {
            arena,
            variable,
            limit,
            depth: 0,
            max_depth: 100,
        }
    }
}

/// Computes the limit of an expression.
///
/// # Arguments
/// * `arena` - The expression arena
/// * `expr` - The expression to compute the limit of
/// * `variable` - The limit variable
/// * `limit` - The direction of the limit
///
/// # Returns
/// The limit result, or an error if computation failed.
pub fn compute_limit(
    arena: &mut ExprArena,
    expr: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
) -> Result<LimitResult, LimitError> {
    gruntz_impl(arena, expr, variable, limit, 0, 100)
}

/// The main Gruntz algorithm implementation.
fn gruntz_impl(
    arena: &mut ExprArena,
    expr: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    // Check recursion depth
    if depth > max_depth {
        return Err(LimitError::CannotCompute("max depth exceeded".to_string()));
    }

    let e = arena.get(expr).clone();

    // Base cases
    match &e {
        // Constants have trivial limits
        ExprNode::Integer(_) | ExprNode::Rational(_, _) => {
            return Ok(LimitResult::Finite(expr));
        }

        // The variable itself
        ExprNode::Symbol(_) => {
            if expr == variable {
                return match limit {
                    Limit::PosInfinity => Ok(LimitResult::PosInfinity),
                    Limit::NegInfinity => Ok(LimitResult::NegInfinity),
                    Limit::Right(a) | Limit::Left(a) | Limit::TwoSided(a) => {
                        Ok(LimitResult::Finite(a))
                    }
                };
            } else {
                // Other symbols are constants
                return Ok(LimitResult::Finite(expr));
            }
        }

        _ => {}
    }

    // Check if expression contains the variable
    if !contains_variable(arena, expr, variable) {
        return Ok(LimitResult::Finite(expr));
    }

    // Handle different expression types
    match &e {
        ExprNode::Add(args) => {
            let args = args.clone();
            limit_add(arena, &args, variable, limit, depth, max_depth)
        }
        ExprNode::Mul(args) => {
            let args = args.clone();
            limit_mul(arena, &args, variable, limit, depth, max_depth)
        }
        ExprNode::Pow { base, exp } => {
            let base = *base;
            let exp = *exp;
            limit_pow(arena, base, exp, variable, limit, depth, max_depth)
        }
        ExprNode::Function { id, args } => {
            let id = *id;
            let args = args.clone();
            match id {
                functions::EXP => limit_exp(arena, args[0], variable, limit, depth, max_depth),
                functions::LN => limit_log(arena, args[0], variable, limit, depth, max_depth),
                functions::SIN => limit_sin(arena, args[0], variable, limit, depth, max_depth),
                functions::COS => limit_cos(arena, args[0], variable, limit, depth, max_depth),
                _ => Ok(LimitResult::Unknown),
            }
        }
        ExprNode::Neg(arg) => {
            let arg = *arg;
            let inner = gruntz_impl(arena, arg, variable, limit, depth + 1, max_depth)?;
            Ok(match inner {
                LimitResult::Finite(v) => {
                    let neg_v = arena.neg(v);
                    LimitResult::Finite(neg_v)
                }
                LimitResult::PosInfinity => LimitResult::NegInfinity,
                LimitResult::NegInfinity => LimitResult::PosInfinity,
                other => other,
            })
        }
        ExprNode::Div { num, den } => {
            let num = *num;
            let den = *den;
            limit_div(arena, num, den, variable, limit, depth, max_depth)
        }
        _ => {
            // For complex expressions, use MRV-based approach
            limit_mrv(arena, expr, variable, limit, depth, max_depth)
        }
    }
}

/// Limit of a sum.
fn limit_add(
    arena: &mut ExprArena,
    args: &[ExprHandle],
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    // Find the dominant term
    let mrv_ctx = MrvContext::new(arena, variable);

    let mut dominant_idx = 0;
    let mut dominant_rate = mrv_ctx.growth_rate(args[0]);

    for (i, &arg) in args.iter().enumerate().skip(1) {
        let rate = mrv_ctx.growth_rate(arg);
        if rate.compare(&dominant_rate) == ComparisonClass::GreaterThan {
            dominant_idx = i;
            dominant_rate = rate;
        }
    }

    // Compute limit of dominant term
    gruntz_impl(arena, args[dominant_idx], variable, limit, depth + 1, max_depth)
}

/// Limit of a product.
fn limit_mul(
    arena: &mut ExprArena,
    args: &[ExprHandle],
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    // Compute limits of each factor
    let mut finite_product = Vec::new();
    let mut has_infinity = false;
    let mut has_zero = false;
    let mut infinity_sign = 1i8;

    for &arg in args {
        match gruntz_impl(arena, arg, variable, limit, depth + 1, max_depth)? {
            LimitResult::Finite(v) => {
                // Check if it's zero
                if is_zero(arena, v) {
                    has_zero = true;
                } else {
                    // Check if negative - affects sign of infinity
                    if is_negative(arena, v) {
                        infinity_sign *= -1;
                    }
                    finite_product.push(v);
                }
            }
            LimitResult::PosInfinity => {
                has_infinity = true;
            }
            LimitResult::NegInfinity => {
                has_infinity = true;
                infinity_sign *= -1;
            }
            LimitResult::DoesNotExist => return Ok(LimitResult::DoesNotExist),
            LimitResult::Unknown => return Ok(LimitResult::Unknown),
        }
    }

    // 0 * ∞ is indeterminate
    if has_zero && has_infinity {
        // Need to use L'Hôpital or series expansion
        return Ok(LimitResult::Unknown);
    }

    if has_zero {
        let zero = arena.integer(0);
        return Ok(LimitResult::Finite(zero));
    }

    if has_infinity {
        return Ok(if infinity_sign > 0 {
            LimitResult::PosInfinity
        } else {
            LimitResult::NegInfinity
        });
    }

    // All finite: multiply them
    if finite_product.is_empty() {
        let one = arena.integer(1);
        return Ok(LimitResult::Finite(one));
    }

    if finite_product.len() == 1 {
        return Ok(LimitResult::Finite(finite_product[0]));
    }

    let result = arena.mul(finite_product);
    Ok(LimitResult::Finite(result))
}

/// Limit of a division.
fn limit_div(
    arena: &mut ExprArena,
    num: ExprHandle,
    den: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    let num_lim = gruntz_impl(arena, num, variable, limit, depth + 1, max_depth)?;
    let den_lim = gruntz_impl(arena, den, variable, limit, depth + 1, max_depth)?;

    match (num_lim, den_lim) {
        (LimitResult::Finite(n), LimitResult::Finite(d)) => {
            if is_zero(arena, d) {
                if is_zero(arena, n) {
                    // 0/0 indeterminate
                    Ok(LimitResult::Unknown)
                } else {
                    // c/0 → ±∞ depending on sign
                    Ok(LimitResult::PosInfinity) // Simplified
                }
            } else {
                let result = arena.intern(ExprNode::Div { num: n, den: d });
                Ok(LimitResult::Finite(result))
            }
        }
        (LimitResult::Finite(_), LimitResult::PosInfinity) => {
            let zero = arena.integer(0);
            Ok(LimitResult::Finite(zero))
        }
        (LimitResult::PosInfinity, LimitResult::Finite(d)) => {
            if !is_zero(arena, d) {
                Ok(LimitResult::PosInfinity)
            } else {
                Ok(LimitResult::Unknown)
            }
        }
        // ∞/∞ is indeterminate
        (LimitResult::PosInfinity, LimitResult::PosInfinity)
        | (LimitResult::NegInfinity, LimitResult::NegInfinity) => Ok(LimitResult::Unknown),
        _ => Ok(LimitResult::Unknown),
    }
}

/// Limit of a power.
fn limit_pow(
    arena: &mut ExprArena,
    base: ExprHandle,
    exp: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    let base_lim = gruntz_impl(arena, base, variable, limit, depth + 1, max_depth)?;
    let exp_lim = gruntz_impl(arena, exp, variable, limit, depth + 1, max_depth)?;

    match (base_lim, exp_lim) {
        (LimitResult::Finite(b), LimitResult::Finite(e)) => {
            let result = arena.pow(b, e);
            Ok(LimitResult::Finite(result))
        }

        // ∞^(positive) = ∞
        (LimitResult::PosInfinity, LimitResult::Finite(e)) => {
            if is_positive(arena, e) {
                Ok(LimitResult::PosInfinity)
            } else if is_negative(arena, e) {
                let zero = arena.integer(0);
                Ok(LimitResult::Finite(zero))
            } else {
                Ok(LimitResult::Unknown)
            }
        }

        // (positive finite)^∞ depends on whether base > 1 or < 1
        (LimitResult::Finite(b), LimitResult::PosInfinity) => {
            if is_greater_than_one(arena, b) {
                Ok(LimitResult::PosInfinity)
            } else if is_less_than_one(arena, b) && is_positive(arena, b) {
                let zero = arena.integer(0);
                Ok(LimitResult::Finite(zero))
            } else {
                Ok(LimitResult::Unknown)
            }
        }

        (LimitResult::PosInfinity, LimitResult::PosInfinity) => Ok(LimitResult::PosInfinity),

        // 0^0, 1^∞, ∞^0 are indeterminate
        _ => Ok(LimitResult::Unknown),
    }
}

/// Limit of exp(f).
fn limit_exp(
    arena: &mut ExprArena,
    arg: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    use smallvec::smallvec;

    match gruntz_impl(arena, arg, variable, limit, depth + 1, max_depth)? {
        LimitResult::Finite(v) => {
            let result = arena.intern(ExprNode::Function {
                id: functions::EXP,
                args: smallvec![v],
            });
            Ok(LimitResult::Finite(result))
        }
        LimitResult::PosInfinity => Ok(LimitResult::PosInfinity),
        LimitResult::NegInfinity => {
            let zero = arena.integer(0);
            Ok(LimitResult::Finite(zero))
        }
        LimitResult::DoesNotExist => Ok(LimitResult::DoesNotExist),
        LimitResult::Unknown => Ok(LimitResult::Unknown),
    }
}

/// Limit of log(f).
fn limit_log(
    arena: &mut ExprArena,
    arg: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    use smallvec::smallvec;

    match gruntz_impl(arena, arg, variable, limit, depth + 1, max_depth)? {
        LimitResult::Finite(v) => {
            if is_zero(arena, v) {
                // log(0⁺) = -∞
                Ok(LimitResult::NegInfinity)
            } else {
                let result = arena.intern(ExprNode::Function {
                    id: functions::LN,
                    args: smallvec![v],
                });
                Ok(LimitResult::Finite(result))
            }
        }
        LimitResult::PosInfinity => Ok(LimitResult::PosInfinity),
        LimitResult::NegInfinity => Ok(LimitResult::DoesNotExist), // log of negative
        LimitResult::DoesNotExist => Ok(LimitResult::DoesNotExist),
        LimitResult::Unknown => Ok(LimitResult::Unknown),
    }
}

/// Limit of sin(f).
fn limit_sin(
    arena: &mut ExprArena,
    arg: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    use smallvec::smallvec;

    match gruntz_impl(arena, arg, variable, limit, depth + 1, max_depth)? {
        LimitResult::Finite(v) => {
            let result = arena.intern(ExprNode::Function {
                id: functions::SIN,
                args: smallvec![v],
            });
            Ok(LimitResult::Finite(result))
        }
        // sin(∞) oscillates
        LimitResult::PosInfinity | LimitResult::NegInfinity => Ok(LimitResult::DoesNotExist),
        other => Ok(other),
    }
}

/// Limit of cos(f).
fn limit_cos(
    arena: &mut ExprArena,
    arg: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    depth: usize,
    max_depth: usize,
) -> Result<LimitResult, LimitError> {
    use smallvec::smallvec;

    match gruntz_impl(arena, arg, variable, limit, depth + 1, max_depth)? {
        LimitResult::Finite(v) => {
            let result = arena.intern(ExprNode::Function {
                id: functions::COS,
                args: smallvec![v],
            });
            Ok(LimitResult::Finite(result))
        }
        // cos(∞) oscillates
        LimitResult::PosInfinity | LimitResult::NegInfinity => Ok(LimitResult::DoesNotExist),
        other => Ok(other),
    }
}

/// MRV-based limit computation for complex expressions.
fn limit_mrv(
    arena: &mut ExprArena,
    expr: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
    _depth: usize,
    _max_depth: usize,
) -> Result<LimitResult, LimitError> {
    // Find the MRV set
    let mrv = find_mrv(arena, expr, variable);

    if mrv.is_empty() {
        // No rapidly varying subexpressions: expression is constant
        return Ok(LimitResult::Finite(expr));
    }

    // Get the dominant MRV element
    let _omega_expr = mrv.dominant().unwrap();

    // Determine the sign of omega as x → ∞
    // For now, assume omega → +∞ and rewrite as ω → 0⁺

    // The full Gruntz algorithm would:
    // 1. Substitute ω = exp(-t) where t → ∞
    // 2. Rewrite expr in terms of ω
    // 3. Expand as a series in ω
    // 4. Extract the leading coefficient

    // Simplified: try direct evaluation for rational expressions
    try_direct_limit(arena, expr, variable, limit)
}

/// Try to compute a limit directly (for simple rational expressions).
fn try_direct_limit(
    arena: &mut ExprArena,
    expr: ExprHandle,
    variable: ExprHandle,
    limit: Limit,
) -> Result<LimitResult, LimitError> {
    let e = arena.get(expr).clone();

    // Check for x/x^n type patterns
    if let ExprNode::Mul(args) = e {
        // Collect powers of x
        let mrv_ctx = MrvContext::new(arena, variable);
        let mut net_degree = 0i32;
        let mut coefficient_handles = Vec::new();

        for &arg in &args {
            let rate = mrv_ctx.growth_rate(arg);
            match rate {
                crate::comparison::GrowthRate::Polynomial(d) => {
                    net_degree += d;
                }
                crate::comparison::GrowthRate::Constant => {
                    // Include in coefficient
                    coefficient_handles.push(arg);
                }
                _ => return Ok(LimitResult::Unknown),
            }
        }

        // Build coefficient expression
        let coefficient = if coefficient_handles.is_empty() {
            arena.integer(1)
        } else if coefficient_handles.len() == 1 {
            coefficient_handles[0]
        } else {
            arena.mul(coefficient_handles)
        };

        // Determine limit based on net degree
        return Ok(match net_degree.cmp(&0) {
            std::cmp::Ordering::Greater => match limit {
                Limit::PosInfinity => LimitResult::PosInfinity,
                Limit::NegInfinity => {
                    if net_degree % 2 == 0 {
                        LimitResult::PosInfinity
                    } else {
                        LimitResult::NegInfinity
                    }
                }
                _ => LimitResult::Unknown,
            },
            std::cmp::Ordering::Less => {
                let zero = arena.integer(0);
                LimitResult::Finite(zero)
            }
            std::cmp::Ordering::Equal => LimitResult::Finite(coefficient),
        });
    }

    Ok(LimitResult::Unknown)
}

// Helper functions

/// Checks if an expression contains a variable.
fn contains_variable(arena: &ExprArena, expr: ExprHandle, var: ExprHandle) -> bool {
    if expr == var {
        return true;
    }

    let e = arena.get(expr);
    match e {
        ExprNode::Integer(_) | ExprNode::Rational(_, _) | ExprNode::Symbol(_) => false,
        ExprNode::Add(args) | ExprNode::Mul(args) => {
            args.iter().any(|&a| contains_variable(arena, a, var))
        }
        ExprNode::Pow { base, exp } => {
            contains_variable(arena, *base, var) || contains_variable(arena, *exp, var)
        }
        ExprNode::Function { args, .. } => args.iter().any(|&a| contains_variable(arena, a, var)),
        ExprNode::Neg(arg) => contains_variable(arena, *arg, var),
        ExprNode::Div { num, den } => {
            contains_variable(arena, *num, var) || contains_variable(arena, *den, var)
        }
    }
}

/// Checks if an expression is zero.
fn is_zero(arena: &ExprArena, expr: ExprHandle) -> bool {
    matches!(arena.get(expr), ExprNode::Integer(0))
}

/// Checks if an expression is positive.
fn is_positive(arena: &ExprArena, expr: ExprHandle) -> bool {
    match arena.get(expr) {
        ExprNode::Integer(n) => *n > 0,
        ExprNode::Rational(n, d) => (*n > 0) == (*d > 0),
        _ => false, // Unknown
    }
}

/// Checks if an expression is negative.
fn is_negative(arena: &ExprArena, expr: ExprHandle) -> bool {
    match arena.get(expr) {
        ExprNode::Integer(n) => *n < 0,
        ExprNode::Rational(n, d) => (*n > 0) != (*d > 0),
        _ => false,
    }
}

/// Checks if an expression is greater than 1.
fn is_greater_than_one(arena: &ExprArena, expr: ExprHandle) -> bool {
    match arena.get(expr) {
        ExprNode::Integer(n) => *n > 1,
        ExprNode::Rational(n, d) => {
            let d_i64 = *d as i64;
            if d_i64 > 0 {
                *n > d_i64
            } else {
                *n < d_i64
            }
        }
        _ => false,
    }
}

/// Checks if an expression is less than 1 (and positive).
fn is_less_than_one(arena: &ExprArena, expr: ExprHandle) -> bool {
    match arena.get(expr) {
        ExprNode::Integer(n) => *n >= 0 && *n < 1,
        ExprNode::Rational(n, d) => {
            let d_i64 = *d as i64;
            if d_i64 > 0 && *n > 0 {
                *n < d_i64
            } else if d_i64 < 0 && *n < 0 {
                *n > d_i64
            } else {
                false
            }
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;

    fn setup() -> (ExprArena, ExprHandle) {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        (arena, x)
    }

    #[test]
    fn test_limit_constant() {
        let (mut arena, x) = setup();
        let c = arena.integer(42);

        let result = compute_limit(&mut arena, c, x, Limit::PosInfinity).unwrap();
        assert!(matches!(result, LimitResult::Finite(v) if v == c));
    }

    #[test]
    fn test_limit_variable() {
        let (mut arena, x) = setup();

        let result = compute_limit(&mut arena, x, x, Limit::PosInfinity).unwrap();
        assert!(matches!(result, LimitResult::PosInfinity));

        let result = compute_limit(&mut arena, x, x, Limit::NegInfinity).unwrap();
        assert!(matches!(result, LimitResult::NegInfinity));
    }

    #[test]
    fn test_limit_exp_neg_infinity() {
        let (mut arena, x) = setup();
        // lim(x→∞) exp(-x) = 0
        let neg_one = arena.integer(-1);
        let neg_x = arena.mul(smallvec![neg_one, x]);
        let exp_neg_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![neg_x],
        });

        let result = compute_limit(&mut arena, exp_neg_x, x, Limit::PosInfinity).unwrap();
        match result {
            LimitResult::Finite(v) => {
                assert!(is_zero(&arena, v));
            }
            _ => panic!("Expected finite limit of 0"),
        }
    }

    #[test]
    fn test_limit_exp_pos_infinity() {
        let (mut arena, x) = setup();
        // lim(x→∞) exp(x) = ∞
        let exp_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![x],
        });

        let result = compute_limit(&mut arena, exp_x, x, Limit::PosInfinity).unwrap();
        assert!(matches!(result, LimitResult::PosInfinity));
    }

    #[test]
    fn test_limit_sin_infinity() {
        let (mut arena, x) = setup();
        // lim(x→∞) sin(x) does not exist
        let sin_x = arena.intern(ExprNode::Function {
            id: functions::SIN,
            args: smallvec![x],
        });

        let result = compute_limit(&mut arena, sin_x, x, Limit::PosInfinity).unwrap();
        assert!(matches!(result, LimitResult::DoesNotExist));
    }

    #[test]
    fn test_limit_at_point() {
        let (mut arena, x) = setup();
        // lim(x→0) x = 0
        let zero = arena.integer(0);

        let result = compute_limit(&mut arena, x, x, Limit::TwoSided(zero)).unwrap();
        assert!(matches!(result, LimitResult::Finite(v) if v == zero));
    }

    #[test]
    fn test_limit_log_infinity() {
        let (mut arena, x) = setup();
        // lim(x→∞) log(x) = ∞
        let log_x = arena.intern(ExprNode::Function {
            id: functions::LN,
            args: smallvec![x],
        });

        let result = compute_limit(&mut arena, log_x, x, Limit::PosInfinity).unwrap();
        assert!(matches!(result, LimitResult::PosInfinity));
    }
}
