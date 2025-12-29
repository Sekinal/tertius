//! Type conversion between expressions and backend types.
//!
//! This module provides bidirectional conversion between the unified
//! `ExprHandle` representation and specialized backend types like
//! `DensePoly<Q>` and `RationalFunction<Q>`.

use smallvec::{smallvec, SmallVec};

use tertius_core::arena::ExprArena;
use tertius_core::expr::{functions, ExprNode, SymbolId};
use tertius_core::handle::ExprHandle;
use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::rationals::Q;
use tertius_rings::traits::Ring;

use crate::rothstein_trager::LogarithmicPart;

/// Bidirectional converter between ExprHandle and backend types.
pub struct ExprConverter<'a> {
    arena: &'a mut ExprArena,
    variable: SymbolId,
    var_handle: ExprHandle,
}

impl<'a> ExprConverter<'a> {
    /// Creates a new converter for the given arena and variable.
    pub fn new(arena: &'a mut ExprArena, variable: SymbolId, var_handle: ExprHandle) -> Self {
        Self {
            arena,
            variable,
            var_handle,
        }
    }

    // === ExprHandle -> Backend Types ===

    /// Converts an expression to a polynomial over Q (if possible).
    pub fn to_polynomial(&self, expr: ExprHandle) -> Option<DensePoly<Q>> {
        self.collect_polynomial(expr)
    }

    /// Converts an expression to a rational function over Q.
    pub fn to_rational_function(&self, expr: ExprHandle) -> Option<RationalFunction<Q>> {
        let node = self.arena.get(expr);

        match node {
            ExprNode::Div { num, den } => {
                let num_poly = self.to_polynomial(*num)?;
                let den_poly = self.to_polynomial(*den)?;
                Some(RationalFunction::new(num_poly, den_poly))
            }
            _ => {
                // Treat as polynomial / 1
                let poly = self.to_polynomial(expr)?;
                Some(RationalFunction::from_poly(poly))
            }
        }
    }

    /// Collects polynomial coefficients from an expression.
    fn collect_polynomial(&self, expr: ExprHandle) -> Option<DensePoly<Q>> {
        let node = self.arena.get(expr);

        match node {
            ExprNode::Integer(n) => Some(DensePoly::constant(Q::from_integer(*n))),

            ExprNode::Rational(num, den) => {
                Some(DensePoly::constant(Q::new(*num, *den as i64)))
            }

            ExprNode::Symbol(id) if *id == self.variable => Some(DensePoly::x()),

            ExprNode::Symbol(_) => {
                // Different variable - treat as constant 0 for now
                // In a full implementation, this would be a multivariate case
                None
            }

            ExprNode::Add(args) => {
                let mut result: DensePoly<Q> = DensePoly::zero();
                for arg in args {
                    let poly = self.collect_polynomial(*arg)?;
                    result = result.add(&poly);
                }
                Some(result)
            }

            ExprNode::Mul(args) => {
                let mut result: DensePoly<Q> = DensePoly::one();
                for arg in args {
                    let poly = self.collect_polynomial(*arg)?;
                    result = result.mul(&poly);
                }
                Some(result)
            }

            ExprNode::Pow { base, exp } => {
                let base_poly = self.collect_polynomial(*base)?;
                let exp_node = self.arena.get(*exp);

                // Only handle non-negative integer exponents
                if let ExprNode::Integer(n) = exp_node {
                    if *n >= 0 {
                        let mut result: DensePoly<Q> = DensePoly::one();
                        for _ in 0..*n {
                            result = result.mul(&base_poly);
                        }
                        return Some(result);
                    }
                }
                None
            }

            ExprNode::Neg(arg) => {
                let poly = self.collect_polynomial(*arg)?;
                Some(poly.neg())
            }

            ExprNode::Div { .. } => None, // Division is not polynomial

            ExprNode::Function { .. } => None, // Functions are not polynomial
        }
    }

    // === Backend Types -> ExprHandle ===

    /// Converts a polynomial back to an ExprHandle.
    pub fn from_polynomial(&mut self, poly: &DensePoly<Q>) -> ExprHandle {
        if poly.is_zero() {
            return self.arena.integer(0);
        }

        let mut result_terms: Vec<ExprHandle> = Vec::new();

        for (i, coeff) in poly.coeffs().iter().enumerate() {
            if !coeff.is_zero() {
                let term = self.make_term(coeff, i);
                result_terms.push(term);
            }
        }

        if result_terms.is_empty() {
            self.arena.integer(0)
        } else if result_terms.len() == 1 {
            result_terms[0]
        } else {
            self.arena.add(SmallVec::from_vec(result_terms))
        }
    }

    /// Converts a rational function back to an ExprHandle.
    pub fn from_rational_function(&mut self, rf: &RationalFunction<Q>) -> ExprHandle {
        let num = self.from_polynomial(rf.numerator());
        let den = self.from_polynomial(rf.denominator());

        // Check if denominator is the constant 1
        if rf.denominator().degree() == 0 && rf.denominator().leading_coeff().is_one() {
            num
        } else {
            self.arena.intern(ExprNode::Div { num, den })
        }
    }

    /// Converts logarithmic terms to an ExprHandle.
    pub fn from_logarithmic_part(&mut self, log_part: &LogarithmicPart<Q>) -> ExprHandle {
        if log_part.is_empty() {
            return self.arena.integer(0);
        }

        let mut result_terms: Vec<ExprHandle> = Vec::new();

        for term in &log_part.terms {
            let coeff_expr = self.from_rational(&term.coefficient);
            let arg_expr = self.from_polynomial(&term.argument);

            // ln(argument)
            let ln_expr = self.arena.intern(ExprNode::Function {
                id: functions::LN,
                args: smallvec![arg_expr],
            });

            // coefficient * ln(argument)
            if term.coefficient.is_one() {
                result_terms.push(ln_expr);
            } else {
                let product = self.arena.mul(smallvec![coeff_expr, ln_expr]);
                result_terms.push(product);
            }
        }

        if result_terms.is_empty() {
            self.arena.integer(0)
        } else if result_terms.len() == 1 {
            result_terms[0]
        } else {
            self.arena.add(SmallVec::from_vec(result_terms))
        }
    }

    /// Converts an integration result with all parts.
    pub fn from_integration_result(
        &mut self,
        poly_part: &DensePoly<Q>,
        rational_part: &RationalFunction<Q>,
        log_part: Option<&LogarithmicPart<Q>>,
    ) -> ExprHandle {
        let mut result_terms: Vec<ExprHandle> = Vec::new();

        // Polynomial part
        if !poly_part.is_zero() {
            result_terms.push(self.from_polynomial(poly_part));
        }

        // Rational part
        if !rational_part.is_zero() {
            result_terms.push(self.from_rational_function(rational_part));
        }

        // Logarithmic part
        if let Some(lp) = log_part {
            if !lp.is_empty() {
                result_terms.push(self.from_logarithmic_part(lp));
            }
        }

        if result_terms.is_empty() {
            self.arena.integer(0)
        } else if result_terms.len() == 1 {
            result_terms[0]
        } else {
            self.arena.add(SmallVec::from_vec(result_terms))
        }
    }

    // === Helper methods ===

    /// Creates a term c * x^n.
    fn make_term(&mut self, coeff: &Q, power: usize) -> ExprHandle {
        let coeff_expr = self.from_rational(coeff);

        if power == 0 {
            return coeff_expr;
        }

        let x_power = if power == 1 {
            self.var_handle
        } else {
            let exp = self.arena.integer(power as i64);
            self.arena.pow(self.var_handle, exp)
        };

        if coeff.is_one() {
            x_power
        } else {
            self.arena.mul(smallvec![coeff_expr, x_power])
        }
    }

    /// Converts a rational number to an ExprHandle.
    fn from_rational(&mut self, q: &Q) -> ExprHandle {
        let (num, den) = q.num_den();

        if den == 1 {
            self.arena.integer(num)
        } else {
            self.arena.intern(ExprNode::Rational(num, den as u64))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_roundtrip() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");

        // Create polynomial: 2x^2 + 3x + 1
        let poly = DensePoly::new(vec![
            Q::from_integer(1),
            Q::from_integer(3),
            Q::from_integer(2),
        ]);

        let mut converter = ExprConverter::new(&mut arena, x_id, x);
        let expr = converter.from_polynomial(&poly);

        // Convert back
        let recovered = converter.to_polynomial(expr);
        assert!(recovered.is_some());
        assert_eq!(recovered.unwrap(), poly);
    }

    #[test]
    fn test_constant_polynomial() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");

        let poly = DensePoly::constant(Q::from_integer(5));

        let mut converter = ExprConverter::new(&mut arena, x_id, x);
        let expr = converter.from_polynomial(&poly);

        // Should just be the integer 5
        let node = arena.get(expr);
        assert!(matches!(node, ExprNode::Integer(5)));
    }

    #[test]
    fn test_to_polynomial_simple() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");
        let two = arena.integer(2);

        // x^2
        let x_squared = arena.pow(x, two);

        let converter = ExprConverter::new(&mut arena, x_id, x);
        let poly = converter.to_polynomial(x_squared);

        assert!(poly.is_some());
        let p = poly.unwrap();
        assert_eq!(p.degree(), 2);
        assert!(p.coeff(0).is_zero());
        assert!(p.coeff(1).is_zero());
        assert!(p.coeff(2).is_one());
    }
}
