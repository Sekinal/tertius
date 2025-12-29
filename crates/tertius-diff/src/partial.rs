//! Partial derivatives and vector calculus operations.
//!
//! This module provides:
//! - Partial derivatives
//! - Gradient (∇f)
//! - Jacobian matrix
//! - Hessian matrix

use tertius_core::arena::ExprArena;
use tertius_core::handle::ExprHandle;

use crate::diff::diff;

/// Computes the partial derivative of an expression with respect to a variable.
///
/// This is equivalent to `diff` but emphasizes that other variables are treated
/// as constants.
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The expression to differentiate
/// * `var` - The variable to differentiate with respect to
///
/// # Returns
///
/// The partial derivative ∂f/∂var.
///
/// # Examples
///
/// ```
/// use tertius_core::arena::ExprArena;
/// use tertius_diff::partial;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
/// let y = arena.symbol("y");
///
/// // f(x, y) = x*y
/// let f = arena.mul(smallvec::smallvec![x, y]);
///
/// // ∂f/∂x = y
/// let df_dx = partial(&mut arena, f, x);
///
/// // ∂f/∂y = x
/// let df_dy = partial(&mut arena, f, y);
/// ```
pub fn partial(arena: &mut ExprArena, expr: ExprHandle, var: ExprHandle) -> ExprHandle {
    diff(arena, expr, var)
}

/// Computes the gradient of a scalar-valued function.
///
/// The gradient is a vector of partial derivatives:
/// ∇f = (∂f/∂x₁, ∂f/∂x₂, ..., ∂f/∂xₙ)
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The scalar expression
/// * `vars` - The variables with respect to which to differentiate
///
/// # Returns
///
/// A vector of partial derivatives.
///
/// # Examples
///
/// ```
/// use tertius_core::arena::ExprArena;
/// use tertius_diff::gradient;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
/// let y = arena.symbol("y");
/// let two = arena.integer(2);
///
/// // f(x, y) = x^2 + y^2
/// let x_sq = arena.pow(x, two);
/// let y_sq = arena.pow(y, two);
/// let f = arena.add(smallvec::smallvec![x_sq, y_sq]);
///
/// // ∇f = (2x, 2y)
/// let grad = gradient(&mut arena, f, &[x, y]);
/// assert_eq!(grad.len(), 2);
/// ```
pub fn gradient(arena: &mut ExprArena, expr: ExprHandle, vars: &[ExprHandle]) -> Vec<ExprHandle> {
    vars.iter().map(|&var| partial(arena, expr, var)).collect()
}

/// Computes the Jacobian matrix of a vector-valued function.
///
/// The Jacobian is a matrix of partial derivatives:
/// J[i,j] = ∂fᵢ/∂xⱼ
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `exprs` - The vector of expressions (f₁, f₂, ..., fₘ)
/// * `vars` - The variables (x₁, x₂, ..., xₙ)
///
/// # Returns
///
/// A matrix (as Vec<Vec<ExprHandle>>) where entry [i][j] is ∂fᵢ/∂xⱼ.
///
/// # Examples
///
/// ```
/// use tertius_core::arena::ExprArena;
/// use tertius_diff::jacobian;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
/// let y = arena.symbol("y");
///
/// // f = (x*y, x+y)
/// let f1 = arena.mul(smallvec::smallvec![x, y]);
/// let f2 = arena.add(smallvec::smallvec![x, y]);
///
/// // J = [[y, x], [1, 1]]
/// let jac = jacobian(&mut arena, &[f1, f2], &[x, y]);
/// assert_eq!(jac.len(), 2);    // 2 rows
/// assert_eq!(jac[0].len(), 2); // 2 columns
/// ```
pub fn jacobian(
    arena: &mut ExprArena,
    exprs: &[ExprHandle],
    vars: &[ExprHandle],
) -> Vec<Vec<ExprHandle>> {
    exprs
        .iter()
        .map(|&expr| vars.iter().map(|&var| partial(arena, expr, var)).collect())
        .collect()
}

/// Computes the Hessian matrix of a scalar-valued function.
///
/// The Hessian is the matrix of second partial derivatives:
/// H[i,j] = ∂²f/∂xᵢ∂xⱼ
///
/// Note: The Hessian is symmetric for smooth functions (Schwarz's theorem).
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The scalar expression
/// * `vars` - The variables (x₁, x₂, ..., xₙ)
///
/// # Returns
///
/// A symmetric matrix where entry [i][j] is ∂²f/∂xᵢ∂xⱼ.
///
/// # Examples
///
/// ```
/// use tertius_core::arena::ExprArena;
/// use tertius_diff::hessian;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
/// let y = arena.symbol("y");
/// let two = arena.integer(2);
///
/// // f(x, y) = x^2 + x*y + y^2
/// let x_sq = arena.pow(x, two);
/// let y_sq = arena.pow(y, two);
/// let xy = arena.mul(smallvec::smallvec![x, y]);
/// let f = arena.add(smallvec::smallvec![x_sq, xy, y_sq]);
///
/// // H = [[2, 1], [1, 2]]
/// let hess = hessian(&mut arena, f, &[x, y]);
/// assert_eq!(hess.len(), 2);
/// ```
pub fn hessian(
    arena: &mut ExprArena,
    expr: ExprHandle,
    vars: &[ExprHandle],
) -> Vec<Vec<ExprHandle>> {
    let n = vars.len();
    let mut result = Vec::with_capacity(n);

    // First, compute the gradient
    let grad = gradient(arena, expr, vars);

    // Then compute the Jacobian of the gradient (which is the Hessian)
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            // H[i,j] = ∂/∂xⱼ (∂f/∂xᵢ) = ∂²f/∂xⱼ∂xᵢ
            row.push(partial(arena, grad[i], vars[j]));
        }
        result.push(row);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_core::expr::ExprNode;

    #[test]
    fn test_partial_derivative() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");

        // f(x, y) = x * y
        let f = arena.mul(smallvec::smallvec![x, y]);

        // ∂f/∂x = y
        let df_dx = partial(&mut arena, f, x);
        match arena.get(df_dx) {
            ExprNode::Symbol(id) => {
                assert_eq!(arena.symbol_name(*id).unwrap(), "y");
            }
            _ => panic!("Expected y"),
        }

        // ∂f/∂y = x
        let df_dy = partial(&mut arena, f, y);
        match arena.get(df_dy) {
            ExprNode::Symbol(id) => {
                assert_eq!(arena.symbol_name(*id).unwrap(), "x");
            }
            _ => panic!("Expected x"),
        }
    }

    #[test]
    fn test_gradient() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let two = arena.integer(2);

        // f(x, y) = x^2 + y^2
        let x_sq = arena.pow(x, two);
        let y_sq = arena.pow(y, two);
        let f = arena.add(smallvec::smallvec![x_sq, y_sq]);

        let grad = gradient(&mut arena, f, &[x, y]);
        assert_eq!(grad.len(), 2);
    }

    #[test]
    fn test_jacobian() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");

        // f = (x*y, x+y)
        let f1 = arena.mul(smallvec::smallvec![x, y]);
        let f2 = arena.add(smallvec::smallvec![x, y]);

        let jac = jacobian(&mut arena, &[f1, f2], &[x, y]);
        assert_eq!(jac.len(), 2);
        assert_eq!(jac[0].len(), 2);
        assert_eq!(jac[1].len(), 2);

        // J[0][0] = ∂(x*y)/∂x = y
        match arena.get(jac[0][0]) {
            ExprNode::Symbol(id) => {
                assert_eq!(arena.symbol_name(*id).unwrap(), "y");
            }
            _ => panic!("Expected y"),
        }

        // J[0][1] = ∂(x*y)/∂y = x
        match arena.get(jac[0][1]) {
            ExprNode::Symbol(id) => {
                assert_eq!(arena.symbol_name(*id).unwrap(), "x");
            }
            _ => panic!("Expected x"),
        }

        // J[1][0] = ∂(x+y)/∂x = 1
        assert!(matches!(arena.get(jac[1][0]), ExprNode::Integer(1)));

        // J[1][1] = ∂(x+y)/∂y = 1
        assert!(matches!(arena.get(jac[1][1]), ExprNode::Integer(1)));
    }

    #[test]
    fn test_hessian() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let two = arena.integer(2);

        // f(x, y) = x^2 + y^2
        let x_sq = arena.pow(x, two);
        let y_sq = arena.pow(y, two);
        let f = arena.add(smallvec::smallvec![x_sq, y_sq]);

        let hess = hessian(&mut arena, f, &[x, y]);
        assert_eq!(hess.len(), 2);
        assert_eq!(hess[0].len(), 2);

        // H[0][0] = ∂²f/∂x² = 2
        assert!(matches!(arena.get(hess[0][0]), ExprNode::Integer(2)));

        // H[0][1] = ∂²f/∂x∂y = 0
        assert!(matches!(arena.get(hess[0][1]), ExprNode::Integer(0)));

        // H[1][0] = ∂²f/∂y∂x = 0 (symmetric)
        assert!(matches!(arena.get(hess[1][0]), ExprNode::Integer(0)));

        // H[1][1] = ∂²f/∂y² = 2
        assert!(matches!(arena.get(hess[1][1]), ExprNode::Integer(2)));
    }
}
