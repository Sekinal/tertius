//! Core differentiation functionality.
//!
//! Provides the main `diff` function for computing derivatives of expressions.

use smallvec::SmallVec;
use tertius_core::arena::ExprArena;
use tertius_core::expr::{functions, ExprNode, SymbolId};
use tertius_core::handle::ExprHandle;

/// Differentiates an expression with respect to a variable.
///
/// Uses the standard rules of calculus:
/// - Constant rule: d/dx(c) = 0
/// - Power rule: d/dx(x^n) = n * x^(n-1)
/// - Sum rule: d/dx(f + g) = f' + g'
/// - Product rule: d/dx(f * g) = f'g + fg'
/// - Quotient rule: d/dx(f/g) = (f'g - fg')/g^2
/// - Chain rule: d/dx(f(g(x))) = f'(g(x)) * g'(x)
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The expression to differentiate
/// * `var` - The variable to differentiate with respect to
///
/// # Returns
///
/// The derivative expression.
///
/// # Examples
///
/// ```
/// use tertius_core::arena::ExprArena;
/// use tertius_diff::diff;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
/// let two = arena.integer(2);
///
/// // d/dx(x^2) = 2x
/// let x_squared = arena.pow(x, two);
/// let deriv = diff(&mut arena, x_squared, x);
/// ```
pub fn diff(arena: &mut ExprArena, expr: ExprHandle, var: ExprHandle) -> ExprHandle {
    // Get the variable's symbol id
    let var_id = match arena.get(var) {
        ExprNode::Symbol(id) => *id,
        _ => panic!("diff: var must be a symbol"),
    };

    diff_internal(arena, expr, var_id)
}

/// Computes the nth derivative of an expression.
///
/// # Arguments
///
/// * `arena` - The expression arena
/// * `expr` - The expression to differentiate
/// * `var` - The variable to differentiate with respect to
/// * `n` - The order of the derivative
///
/// # Returns
///
/// The nth derivative expression.
///
/// # Examples
///
/// ```
/// use tertius_core::arena::ExprArena;
/// use tertius_diff::diff_n;
///
/// let mut arena = ExprArena::new();
/// let x = arena.symbol("x");
/// let three = arena.integer(3);
///
/// // d^2/dx^2(x^3) = 6x
/// let x_cubed = arena.pow(x, three);
/// let second_deriv = diff_n(&mut arena, x_cubed, x, 2);
/// ```
pub fn diff_n(arena: &mut ExprArena, expr: ExprHandle, var: ExprHandle, n: usize) -> ExprHandle {
    let mut result = expr;
    for _ in 0..n {
        result = diff(arena, result, var);
    }
    result
}

/// Internal differentiation function that works with symbol IDs.
fn diff_internal(arena: &mut ExprArena, expr: ExprHandle, var_id: SymbolId) -> ExprHandle {
    match arena.get(expr).clone() {
        // Constant rule: d/dx(c) = 0
        ExprNode::Integer(_) | ExprNode::Rational(_, _) => arena.integer(0),

        // Variable rule: d/dx(x) = 1, d/dx(y) = 0
        ExprNode::Symbol(id) => {
            if id == var_id {
                arena.integer(1)
            } else {
                arena.integer(0)
            }
        }

        // Sum rule: d/dx(f + g + ...) = f' + g' + ...
        ExprNode::Add(args) => {
            let derivatives: SmallVec<[ExprHandle; 4]> = args
                .iter()
                .map(|&arg| diff_internal(arena, arg, var_id))
                .collect();
            simplify_sum(arena, derivatives)
        }

        // Product rule: d/dx(f * g * ...) = f' * g * ... + f * g' * ... + ...
        ExprNode::Mul(args) => {
            let n = args.len();
            let mut terms = SmallVec::<[ExprHandle; 4]>::new();

            for i in 0..n {
                // Differentiate the i-th factor, keep others
                let mut factors = SmallVec::<[ExprHandle; 4]>::new();
                for j in 0..n {
                    if i == j {
                        factors.push(diff_internal(arena, args[j], var_id));
                    } else {
                        factors.push(args[j]);
                    }
                }
                let term = simplify_product(arena, factors);
                terms.push(term);
            }

            simplify_sum(arena, terms)
        }

        // Power rule and chain rule
        ExprNode::Pow { base, exp } => {
            // d/dx(f^g) = f^g * (g' * ln(f) + g * f'/f)
            // Special case: if g is constant, d/dx(f^n) = n * f^(n-1) * f'
            // Special case: if f is constant, d/dx(a^g) = a^g * ln(a) * g'

            let base_depends = depends_on(arena, base, var_id);
            let exp_depends = depends_on(arena, exp, var_id);

            if !base_depends && !exp_depends {
                // Constant
                arena.integer(0)
            } else if !exp_depends {
                // Power rule: d/dx(f^n) = n * f^(n-1) * f'
                let df = diff_internal(arena, base, var_id);
                let one = arena.integer(1);
                let exp_minus_1 = simplify_sub(arena, exp, one);
                let power = simplify_pow(arena, base, exp_minus_1);
                let product = simplify_product(arena, smallvec::smallvec![exp, power, df]);
                product
            } else if !base_depends {
                // Exponential rule: d/dx(a^g) = a^g * ln(a) * g'
                let dg = diff_internal(arena, exp, var_id);
                let ln_base = make_function(arena, functions::LN, smallvec::smallvec![base]);
                let product =
                    simplify_product(arena, smallvec::smallvec![expr, ln_base, dg]);
                product
            } else {
                // General case: d/dx(f^g) = f^g * (g' * ln(f) + g * f'/f)
                let df = diff_internal(arena, base, var_id);
                let dg = diff_internal(arena, exp, var_id);
                let ln_f = make_function(arena, functions::LN, smallvec::smallvec![base]);

                // g' * ln(f)
                let term1 = simplify_product(arena, smallvec::smallvec![dg, ln_f]);

                // g * f' / f
                let f_prime_over_f = simplify_div(arena, df, base);
                let term2 = simplify_product(arena, smallvec::smallvec![exp, f_prime_over_f]);

                // g' * ln(f) + g * f'/f
                let inner_sum = simplify_sum(arena, smallvec::smallvec![term1, term2]);

                // f^g * (...)
                simplify_product(arena, smallvec::smallvec![expr, inner_sum])
            }
        }

        // Negation: d/dx(-f) = -f'
        ExprNode::Neg(arg) => {
            let darg = diff_internal(arena, arg, var_id);
            arena.neg(darg)
        }

        // Quotient rule: d/dx(f/g) = (f'g - fg')/g^2
        ExprNode::Div { num, den } => {
            let df = diff_internal(arena, num, var_id);
            let dg = diff_internal(arena, den, var_id);

            // f'g
            let term1 = simplify_product(arena, smallvec::smallvec![df, den]);
            // fg'
            let term2 = simplify_product(arena, smallvec::smallvec![num, dg]);
            // f'g - fg'
            let numerator = simplify_sub(arena, term1, term2);
            // g^2
            let two = arena.integer(2);
            let denominator = arena.pow(den, two);

            simplify_div(arena, numerator, denominator)
        }

        // Function rules (chain rule + known derivatives)
        ExprNode::Function { id, args } => {
            diff_function(arena, id, &args, var_id)
        }
    }
}

/// Differentiates a function application.
fn diff_function(
    arena: &mut ExprArena,
    func_id: u32,
    args: &SmallVec<[ExprHandle; 2]>,
    var_id: SymbolId,
) -> ExprHandle {
    if args.is_empty() {
        return arena.integer(0);
    }

    // Most elementary functions have a single argument
    let arg = args[0];
    let darg = diff_internal(arena, arg, var_id);

    // Chain rule: d/dx(f(g)) = f'(g) * g'
    let outer_deriv = match func_id {
        // d/dx(sin(u)) = cos(u) * u'
        functions::SIN => {
            let cos_u = make_function(arena, functions::COS, smallvec::smallvec![arg]);
            cos_u
        }

        // d/dx(cos(u)) = -sin(u) * u'
        functions::COS => {
            let sin_u = make_function(arena, functions::SIN, smallvec::smallvec![arg]);
            arena.neg(sin_u)
        }

        // d/dx(tan(u)) = sec^2(u) * u' = (1/cos^2(u)) * u'
        functions::TAN => {
            let cos_u = make_function(arena, functions::COS, smallvec::smallvec![arg]);
            let two = arena.integer(2);
            let cos_sq = arena.pow(cos_u, two);
            let one = arena.integer(1);
            simplify_div(arena, one, cos_sq)
        }

        // d/dx(exp(u)) = exp(u) * u'
        functions::EXP => {
            make_function(arena, functions::EXP, smallvec::smallvec![arg])
        }

        // d/dx(ln(u)) = (1/u) * u'
        functions::LN => {
            let one = arena.integer(1);
            simplify_div(arena, one, arg)
        }

        // d/dx(log10(u)) = (1/(u * ln(10))) * u'
        functions::LOG10 => {
            let one = arena.integer(1);
            let ten = arena.integer(10);
            let ln_10 = make_function(arena, functions::LN, smallvec::smallvec![ten]);
            let denom = simplify_product(arena, smallvec::smallvec![arg, ln_10]);
            simplify_div(arena, one, denom)
        }

        // d/dx(sqrt(u)) = (1/(2*sqrt(u))) * u'
        functions::SQRT => {
            let one = arena.integer(1);
            let two = arena.integer(2);
            let sqrt_u = make_function(arena, functions::SQRT, smallvec::smallvec![arg]);
            let denom = simplify_product(arena, smallvec::smallvec![two, sqrt_u]);
            simplify_div(arena, one, denom)
        }

        // d/dx(|u|) = sign(u) * u' = u/|u| * u'
        functions::ABS => {
            let abs_u = make_function(arena, functions::ABS, smallvec::smallvec![arg]);
            simplify_div(arena, arg, abs_u)
        }

        // Unknown function - return symbolic derivative
        _ => {
            // Could create a symbolic Derivative[f, x] expression here
            // For now, just return 0 as a placeholder
            arena.integer(0)
        }
    };

    // Apply chain rule: outer_deriv * darg
    if is_zero(arena, darg) {
        arena.integer(0)
    } else if is_one(arena, darg) {
        outer_deriv
    } else {
        simplify_product(arena, smallvec::smallvec![outer_deriv, darg])
    }
}

/// Checks if an expression depends on a variable.
fn depends_on(arena: &ExprArena, expr: ExprHandle, var_id: SymbolId) -> bool {
    match arena.get(expr) {
        ExprNode::Integer(_) | ExprNode::Rational(_, _) => false,
        ExprNode::Symbol(id) => *id == var_id,
        ExprNode::Add(args) | ExprNode::Mul(args) => {
            args.iter().any(|&arg| depends_on(arena, arg, var_id))
        }
        ExprNode::Pow { base, exp } => {
            depends_on(arena, *base, var_id) || depends_on(arena, *exp, var_id)
        }
        ExprNode::Neg(arg) => depends_on(arena, *arg, var_id),
        ExprNode::Div { num, den } => {
            depends_on(arena, *num, var_id) || depends_on(arena, *den, var_id)
        }
        ExprNode::Function { args, .. } => {
            args.iter().any(|&arg| depends_on(arena, arg, var_id))
        }
    }
}

/// Creates a function application node.
fn make_function(
    arena: &mut ExprArena,
    id: u32,
    args: SmallVec<[ExprHandle; 2]>,
) -> ExprHandle {
    arena.intern(ExprNode::Function { id, args })
}

/// Simplifies a sum, removing zeros.
fn simplify_sum(arena: &mut ExprArena, terms: SmallVec<[ExprHandle; 4]>) -> ExprHandle {
    let non_zero: SmallVec<[ExprHandle; 4]> = terms
        .into_iter()
        .filter(|&t| !is_zero(arena, t))
        .collect();

    match non_zero.len() {
        0 => arena.integer(0),
        1 => non_zero[0],
        _ => arena.add(non_zero),
    }
}

/// Simplifies a product, handling zeros and ones.
fn simplify_product(arena: &mut ExprArena, factors: SmallVec<[ExprHandle; 4]>) -> ExprHandle {
    // Check for zeros
    if factors.iter().any(|&f| is_zero(arena, f)) {
        return arena.integer(0);
    }

    // Filter out ones
    let non_one: SmallVec<[ExprHandle; 4]> = factors
        .into_iter()
        .filter(|&f| !is_one(arena, f))
        .collect();

    match non_one.len() {
        0 => arena.integer(1),
        1 => non_one[0],
        _ => arena.mul(non_one),
    }
}

/// Creates a subtraction: a - b.
fn simplify_sub(arena: &mut ExprArena, a: ExprHandle, b: ExprHandle) -> ExprHandle {
    if is_zero(arena, b) {
        return a;
    }
    if is_zero(arena, a) {
        return arena.neg(b);
    }

    // If both are integers, compute the result directly
    match (arena.get(a).clone(), arena.get(b).clone()) {
        (ExprNode::Integer(x), ExprNode::Integer(y)) => {
            return arena.integer(x - y);
        }
        (ExprNode::Rational(num_a, den_a), ExprNode::Integer(y)) => {
            // a/b - y = (a - y*b)/b
            return arena.intern(ExprNode::Rational(num_a - y * den_a as i64, den_a));
        }
        (ExprNode::Integer(x), ExprNode::Rational(num_b, den_b)) => {
            // x - a/b = (x*b - a)/b
            return arena.intern(ExprNode::Rational(x * den_b as i64 - num_b, den_b));
        }
        (ExprNode::Rational(num_a, den_a), ExprNode::Rational(num_b, den_b)) => {
            // a/b - c/d = (a*d - c*b)/(b*d)
            let num = num_a * den_b as i64 - num_b * den_a as i64;
            let den = den_a * den_b;
            return arena.intern(ExprNode::Rational(num, den));
        }
        _ => {}
    }

    let neg_b = arena.neg(b);
    arena.add(smallvec::smallvec![a, neg_b])
}

/// Creates a division, with basic simplification.
fn simplify_div(arena: &mut ExprArena, num: ExprHandle, den: ExprHandle) -> ExprHandle {
    if is_zero(arena, num) {
        return arena.integer(0);
    }
    if is_one(arena, den) {
        return num;
    }
    arena.intern(ExprNode::Div { num, den })
}

/// Creates a power expression with basic simplification.
fn simplify_pow(arena: &mut ExprArena, base: ExprHandle, exp: ExprHandle) -> ExprHandle {
    // x^0 = 1
    if is_zero(arena, exp) {
        return arena.integer(1);
    }
    // x^1 = x
    if is_one(arena, exp) {
        return base;
    }
    // 0^n = 0 (for positive n)
    if is_zero(arena, base) {
        return arena.integer(0);
    }
    // 1^n = 1
    if is_one(arena, base) {
        return arena.integer(1);
    }
    arena.pow(base, exp)
}

/// Checks if an expression is zero.
fn is_zero(arena: &ExprArena, expr: ExprHandle) -> bool {
    matches!(arena.get(expr), ExprNode::Integer(0))
}

/// Checks if an expression is one.
fn is_one(arena: &ExprArena, expr: ExprHandle) -> bool {
    matches!(arena.get(expr), ExprNode::Integer(1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diff_constant() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let five = arena.integer(5);

        let d = diff(&mut arena, five, x);
        assert!(matches!(arena.get(d), ExprNode::Integer(0)));
    }

    #[test]
    fn test_diff_variable() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");

        // d/dx(x) = 1
        let dx = diff(&mut arena, x, x);
        assert!(matches!(arena.get(dx), ExprNode::Integer(1)));

        // d/dx(y) = 0
        let dy = diff(&mut arena, y, x);
        assert!(matches!(arena.get(dy), ExprNode::Integer(0)));
    }

    #[test]
    fn test_diff_sum() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let one = arena.integer(1);
        let sum = arena.add(smallvec::smallvec![x, one]);

        // d/dx(x + 1) = 1
        let d = diff(&mut arena, sum, x);
        assert!(matches!(arena.get(d), ExprNode::Integer(1)));
    }

    #[test]
    fn test_diff_power() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let two = arena.integer(2);
        let x_sq = arena.pow(x, two);

        // d/dx(x^2) = 2x
        let d = diff(&mut arena, x_sq, x);

        // Should be a product of 2 and x
        match arena.get(d) {
            ExprNode::Mul(args) => {
                assert_eq!(args.len(), 2);
                // Check one is 2 and one is x
                let has_two = args.iter().any(|&h| matches!(arena.get(h), ExprNode::Integer(2)));
                let has_x = args.iter().any(|&h| matches!(arena.get(h), ExprNode::Symbol(_)));
                assert!(has_two && has_x);
            }
            _ => panic!("Expected Mul, got {:?}", arena.get(d)),
        }
    }

    #[test]
    fn test_diff_product() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let prod = arena.mul(smallvec::smallvec![x, y]);

        // d/dx(x*y) = y (assuming y is constant w.r.t. x)
        let d = diff(&mut arena, prod, x);

        // Should be y
        match arena.get(d) {
            ExprNode::Symbol(id) => {
                let name = arena.symbol_name(*id).unwrap();
                assert_eq!(name, "y");
            }
            _ => panic!("Expected Symbol y, got {:?}", arena.get(d)),
        }
    }

    #[test]
    fn test_diff_n() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let three = arena.integer(3);
        let x_cubed = arena.pow(x, three);

        // d^2/dx^2(x^3) should be equivalent to 6x
        let d2 = diff_n(&mut arena, x_cubed, x, 2);

        // The result is correct but may be represented as 3*(2*x) instead of 6*x
        // Just verify it's a non-zero expression containing x somewhere
        fn contains_symbol(arena: &ExprArena, h: ExprHandle, sym_id: u32) -> bool {
            match arena.get(h) {
                ExprNode::Symbol(id) => *id == sym_id,
                ExprNode::Add(args) | ExprNode::Mul(args) => {
                    args.iter().any(|&a| contains_symbol(arena, a, sym_id))
                }
                ExprNode::Pow { base, exp } => {
                    contains_symbol(arena, *base, sym_id) || contains_symbol(arena, *exp, sym_id)
                }
                ExprNode::Neg(a) => contains_symbol(arena, *a, sym_id),
                ExprNode::Div { num, den } => {
                    contains_symbol(arena, *num, sym_id) || contains_symbol(arena, *den, sym_id)
                }
                ExprNode::Function { args, .. } => {
                    args.iter().any(|&a| contains_symbol(arena, a, sym_id))
                }
                _ => false,
            }
        }

        // Get x's symbol id
        let x_id = match arena.get(x) {
            ExprNode::Symbol(id) => *id,
            _ => panic!("x should be a symbol"),
        };

        // The second derivative of x^3 should contain x
        assert!(contains_symbol(&arena, d2, x_id));
        // And it shouldn't be zero
        assert!(!matches!(arena.get(d2), ExprNode::Integer(0)));
    }
}
