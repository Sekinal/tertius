//! Public API for unified integration.
//!
//! This module provides the main entry points for integration that work
//! like Mathematica's `Integrate[]` function.
//!
//! # Example
//!
//! ```ignore
//! use tertius_integrate::integrate;
//! use tertius_core::ExprArena;
//!
//! let mut arena = ExprArena::new();
//! let x = arena.symbol("x");
//! let x_squared = arena.pow(x, arena.integer(2));
//!
//! let result = integrate(&mut arena, x_squared, x);
//! // result contains x^3/3
//! ```

use tertius_core::arena::ExprArena;
use tertius_core::handle::ExprHandle;

use super::dispatch::{IntegrationDispatcher, IntegrationOptions};
use super::result::{IntegrationResult, NumericalResult};
use crate::numerical::adaptive_integrate;

/// Computes the indefinite integral of an expression.
///
/// This is the main entry point for symbolic integration. It automatically
/// classifies the integrand and dispatches to the appropriate algorithm.
///
/// # Arguments
///
/// * `arena` - The expression arena (will be mutated to add result nodes)
/// * `expr` - The integrand expression
/// * `var` - The integration variable (must be a Symbol)
///
/// # Returns
///
/// An `IntegrationResult` containing:
/// - `Symbolic`: A closed-form antiderivative
/// - `SpecialFunction`: Result in terms of special functions
/// - `Elliptic`: Elliptic integral form
/// - `NonElementary`: Proof of non-integrability
/// - `Unknown`: Could not determine
///
/// # Example
///
/// ```ignore
/// use tertius_integrate::integrate;
/// use tertius_core::ExprArena;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
///
/// // ∫ x² dx = x³/3
/// let x_squared = arena.pow(x, arena.integer(2));
/// let result = integrate(&mut arena, x_squared, x);
///
/// // ∫ sin(x) dx = -cos(x)
/// let sin_x = arena.intern(ExprNode::Function {
///     id: functions::SIN,
///     args: smallvec![x],
/// });
/// let result = integrate(&mut arena, sin_x, x);
/// ```
pub fn integrate(arena: &mut ExprArena, expr: ExprHandle, var: ExprHandle) -> IntegrationResult {
    let mut dispatcher = IntegrationDispatcher::new(arena);
    dispatcher.integrate(expr, var)
}

/// Computes the indefinite integral with custom options.
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The integrand expression
/// * `var` - The integration variable
/// * `options` - Custom integration options
///
/// # Example
///
/// ```ignore
/// use tertius_integrate::{integrate_with_options, IntegrationOptions};
///
/// let options = IntegrationOptions {
///     verify: true,
///     heuristics_first: false,
///     ..Default::default()
/// };
///
/// let result = integrate_with_options(&mut arena, expr, x, options);
/// ```
pub fn integrate_with_options(
    arena: &mut ExprArena,
    expr: ExprHandle,
    var: ExprHandle,
    options: IntegrationOptions,
) -> IntegrationResult {
    let mut dispatcher = IntegrationDispatcher::with_options(arena, options);
    dispatcher.integrate(expr, var)
}

/// Computes a definite integral.
///
/// First attempts symbolic integration, then evaluates F(upper) - F(lower).
/// Falls back to numerical integration if symbolic methods fail.
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The integrand
/// * `var` - The integration variable
/// * `lower` - Lower bound of integration
/// * `upper` - Upper bound of integration
///
/// # Example
///
/// ```ignore
/// use tertius_integrate::integrate_definite;
///
/// // ∫₀¹ x² dx = 1/3
/// let result = integrate_definite(
///     &mut arena,
///     x_squared,
///     x,
///     arena.integer(0),
///     arena.integer(1),
/// );
/// ```
pub fn integrate_definite(
    arena: &mut ExprArena,
    expr: ExprHandle,
    var: ExprHandle,
    lower: ExprHandle,
    upper: ExprHandle,
) -> IntegrationResult {
    let mut dispatcher = IntegrationDispatcher::new(arena);
    dispatcher.integrate_definite(expr, var, lower, upper)
}

/// Computes a definite integral with custom options.
pub fn integrate_definite_with_options(
    arena: &mut ExprArena,
    expr: ExprHandle,
    var: ExprHandle,
    lower: ExprHandle,
    upper: ExprHandle,
    options: IntegrationOptions,
) -> IntegrationResult {
    let mut dispatcher = IntegrationDispatcher::with_options(arena, options);
    dispatcher.integrate_definite(expr, var, lower, upper)
}

/// Computes a purely numerical definite integral.
///
/// This bypasses symbolic methods entirely and uses adaptive Gauss-Kronrod
/// quadrature for direct numerical integration.
///
/// # Arguments
///
/// * `f` - A closure that evaluates the integrand at a point
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `tolerance` - Desired absolute and relative tolerance
///
/// # Returns
///
/// A `NumericalResult` containing the value and error estimate.
///
/// # Example
///
/// ```ignore
/// use tertius_integrate::integrate_numerical;
///
/// // ∫₀¹ x² dx
/// let result = integrate_numerical(&|x| x * x, 0.0, 1.0, 1e-10);
/// assert!((result.value - 1.0/3.0).abs() < 1e-8);
/// ```
pub fn integrate_numerical<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    tolerance: f64,
) -> NumericalResult {
    let result = adaptive_integrate(f, a, b, tolerance, tolerance, 1000);

    NumericalResult::new(result.value, result.error)
        .with_evaluations(result.evaluations)
        .with_converged(result.converged)
}

/// Computes a numerical integral with custom parameters.
///
/// # Arguments
///
/// * `f` - The integrand function
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `abs_tol` - Absolute tolerance
/// * `rel_tol` - Relative tolerance
/// * `max_subdivisions` - Maximum number of interval subdivisions
pub fn integrate_numerical_with_params<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    abs_tol: f64,
    rel_tol: f64,
    max_subdivisions: usize,
) -> NumericalResult {
    let result = adaptive_integrate(f, a, b, abs_tol, rel_tol, max_subdivisions);

    NumericalResult::new(result.value, result.error)
        .with_evaluations(result.evaluations)
        .with_converged(result.converged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_core::expr::{functions, ExprNode};
    use smallvec::smallvec;

    #[test]
    fn test_integrate_x_squared() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let two = arena.integer(2);
        let x_squared = arena.pow(x, two);

        let result = integrate(&mut arena, x_squared, x);
        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_constant() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let five = arena.integer(5);

        let result = integrate(&mut arena, five, x);
        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_sin() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let sin_x = arena.intern(ExprNode::Function {
            id: functions::SIN,
            args: smallvec![x],
        });

        let result = integrate(&mut arena, sin_x, x);
        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_cos() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let cos_x = arena.intern(ExprNode::Function {
            id: functions::COS,
            args: smallvec![x],
        });

        let result = integrate(&mut arena, cos_x, x);
        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_exp() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let exp_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![x],
        });

        let result = integrate(&mut arena, exp_x, x);
        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_one_over_x() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let one = arena.integer(1);
        let one_over_x = arena.intern(ExprNode::Div { num: one, den: x });

        let result = integrate(&mut arena, one_over_x, x);
        assert!(result.is_symbolic());
    }

    #[test]
    fn test_integrate_numerical() {
        // ∫₀¹ x² dx = 1/3
        let result = integrate_numerical(&|x| x * x, 0.0, 1.0, 1e-10);
        assert!((result.value - 1.0 / 3.0).abs() < 1e-8);
        assert!(result.converged);
    }

    #[test]
    fn test_integrate_numerical_sin() {
        // ∫₀^π sin(x) dx = 2
        let result = integrate_numerical(&|x| x.sin(), 0.0, std::f64::consts::PI, 1e-10);
        assert!((result.value - 2.0).abs() < 1e-8);
    }
}
