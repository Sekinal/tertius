//! Symbolic differentiation for Tertius CAS.
//!
//! This crate provides symbolic differentiation of expressions,
//! including partial derivatives, gradients, Jacobians, and Hessians.
//!
//! # Examples
//!
//! ```
//! use tertius_core::arena::ExprArena;
//! use tertius_diff::{diff, diff_n};
//!
//! let mut arena = ExprArena::new();
//! let x = arena.symbol("x");
//!
//! // Create x^2
//! let two = arena.integer(2);
//! let x_squared = arena.pow(x, two);
//!
//! // Differentiate: d/dx (x^2) = 2x
//! let derivative = diff(&mut arena, x_squared, x);
//! ```

mod diff;
mod partial;
mod rules;
mod vector;

pub use diff::{diff, diff_n, diff_with_assumptions};
pub use partial::{gradient, hessian, jacobian, partial};
pub use vector::{
    curl_cartesian, curl_cylindrical, curl_spherical, divergence_cartesian, divergence_cylindrical,
    divergence_spherical, gradient_cartesian, gradient_cylindrical, gradient_spherical,
    laplacian_cartesian, laplacian_cylindrical, laplacian_spherical,
};

#[cfg(test)]
mod tests {
    use tertius_core::assumptions::{Assumption, AssumptionSet};
    use tertius_core::expr::ExprNode;
    use tertius_core::ExprArena;

    use crate::{diff, diff_with_assumptions};

    #[test]
    fn test_diff_with_assumptions_matches_diff() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let two = arena.integer(2);
        let x_sq = arena.pow(x, two);

        let mut assumptions = AssumptionSet::new();
        assumptions.assume("x", Assumption::Real);

        let plain = diff(&mut arena, x_sq, x);
        let with_assumptions = diff_with_assumptions(&mut arena, x_sq, x, &assumptions);

        assert_eq!(plain, with_assumptions);
        assert!(!matches!(arena.get(with_assumptions), ExprNode::Integer(0)));
    }
}
