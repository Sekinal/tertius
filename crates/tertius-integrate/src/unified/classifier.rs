//! Expression classifier for integration.
//!
//! Analyzes expressions to determine their structure and select the
//! appropriate integration algorithm.

use tertius_core::arena::ExprArena;
use tertius_core::expr::{functions, ExprNode, FunctionId, SymbolId};
use tertius_core::handle::ExprHandle;

/// Classification of an integrand's structure.
#[derive(Clone, Debug, PartialEq)]
pub enum IntegrandClass {
    /// Constant expression (no variable dependence).
    Constant,

    /// Polynomial in the integration variable.
    Polynomial {
        /// Degree of the polynomial.
        degree: usize,
    },

    /// Rational function P(x)/Q(x).
    Rational {
        /// Degree of numerator.
        num_degree: usize,
        /// Degree of denominator.
        den_degree: usize,
        /// Whether algebraic number field extension may be needed.
        needs_algebraic: bool,
    },

    /// Algebraic function involving radicals.
    Algebraic {
        /// Genus of the associated curve.
        genus: AlgebraicGenus,
        /// Degree of the radicand.
        radicand_degree: usize,
    },

    /// Transcendental involving log/exp/trig.
    Transcendental {
        /// Profile of transcendental functions present.
        profile: TranscendentalProfile,
    },

    /// Mixed algebraic and transcendental.
    Mixed {
        /// Algebraic part.
        has_algebraic: bool,
        /// Transcendental profile.
        profile: TranscendentalProfile,
    },

    /// Unclassified or too complex.
    Unknown,
}

/// Genus of an algebraic curve, determines integration method.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlgebraicGenus {
    /// Genus 0: can be rationalized (deg ≤ 2)
    Rational,
    /// Genus 1: elliptic integrals (deg 3-4)
    Elliptic,
    /// Genus ≥ 2: hyperelliptic, non-elementary
    Hyperelliptic(usize),
}

impl AlgebraicGenus {
    /// Computes the genus from radicand degree.
    pub fn from_degree(deg: usize) -> Self {
        match deg {
            0..=2 => AlgebraicGenus::Rational,
            3 | 4 => AlgebraicGenus::Elliptic,
            n => AlgebraicGenus::Hyperelliptic((n - 1) / 2),
        }
    }
}

/// Profile of transcendental functions in an expression.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TranscendentalProfile {
    /// Contains exp or e^x terms.
    pub has_exp: bool,
    /// Contains ln or log terms.
    pub has_log: bool,
    /// Contains sin, cos, tan, etc.
    pub has_trig: bool,
    /// Contains arcsin, arccos, arctan, etc.
    pub has_inverse_trig: bool,
    /// Maximum nesting depth of transcendental functions.
    pub nesting_depth: usize,
    /// List of function IDs present.
    pub functions: Vec<FunctionId>,
}

impl TranscendentalProfile {
    /// Returns true if any transcendental functions are present.
    pub fn has_any(&self) -> bool {
        self.has_exp || self.has_log || self.has_trig || self.has_inverse_trig
    }

    /// Returns true if only logarithms are present (Risch-friendly).
    pub fn is_logarithmic_only(&self) -> bool {
        self.has_log && !self.has_exp && !self.has_trig && !self.has_inverse_trig
    }

    /// Returns true if only exponentials are present (Risch-friendly).
    pub fn is_exponential_only(&self) -> bool {
        self.has_exp && !self.has_log && !self.has_trig && !self.has_inverse_trig
    }
}

/// Analyzes an expression to determine its integrand class.
pub struct IntegrandAnalyzer<'a> {
    arena: &'a ExprArena,
    variable: SymbolId,
}

impl<'a> IntegrandAnalyzer<'a> {
    /// Creates a new analyzer for the given arena and integration variable.
    pub fn new(arena: &'a ExprArena, variable: SymbolId) -> Self {
        Self { arena, variable }
    }

    /// Classifies the integrand.
    pub fn classify(&self, expr: ExprHandle) -> IntegrandClass {
        // Check if constant first
        if !self.depends_on_variable(expr) {
            return IntegrandClass::Constant;
        }

        // Analyze structure
        self.analyze_structure(expr)
    }

    /// Checks if expression depends on the integration variable.
    pub fn depends_on_variable(&self, expr: ExprHandle) -> bool {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Symbol(id) => *id == self.variable,
            ExprNode::Integer(_) | ExprNode::Rational(_, _) => false,
            ExprNode::Add(args) | ExprNode::Mul(args) => {
                args.iter().any(|h| self.depends_on_variable(*h))
            }
            ExprNode::Pow { base, exp } => {
                self.depends_on_variable(*base) || self.depends_on_variable(*exp)
            }
            ExprNode::Neg(arg) => self.depends_on_variable(*arg),
            ExprNode::Div { num, den } => {
                self.depends_on_variable(*num) || self.depends_on_variable(*den)
            }
            ExprNode::Function { args, .. } => {
                args.iter().any(|h| self.depends_on_variable(*h))
            }
        }
    }

    /// Analyzes the structure of an expression.
    fn analyze_structure(&self, expr: ExprHandle) -> IntegrandClass {
        let trans_profile = self.transcendental_profile(expr, 0);
        let has_transcendental = trans_profile.has_any();
        let has_radical = self.has_radical(expr);

        if !has_transcendental && !has_radical {
            // Pure polynomial or rational
            if let Some((num_deg, den_deg)) = self.rational_degrees(expr) {
                if den_deg == 0 {
                    IntegrandClass::Polynomial { degree: num_deg }
                } else {
                    IntegrandClass::Rational {
                        num_degree: num_deg,
                        den_degree: den_deg,
                        needs_algebraic: den_deg > 2,
                    }
                }
            } else {
                IntegrandClass::Unknown
            }
        } else if has_radical && !has_transcendental {
            // Algebraic function
            let rad_deg = self.radicand_degree(expr);
            let genus = AlgebraicGenus::from_degree(rad_deg);
            IntegrandClass::Algebraic {
                genus,
                radicand_degree: rad_deg,
            }
        } else if has_transcendental && !has_radical {
            // Pure transcendental
            IntegrandClass::Transcendental {
                profile: trans_profile,
            }
        } else {
            // Mixed algebraic-transcendental
            IntegrandClass::Mixed {
                has_algebraic: has_radical,
                profile: trans_profile,
            }
        }
    }

    /// Builds a profile of transcendental functions in the expression.
    fn transcendental_profile(&self, expr: ExprHandle, depth: usize) -> TranscendentalProfile {
        let node = self.arena.get(expr);
        let mut profile = TranscendentalProfile::default();
        profile.nesting_depth = depth;

        match node {
            ExprNode::Integer(_) | ExprNode::Rational(_, _) | ExprNode::Symbol(_) => {}

            ExprNode::Add(args) | ExprNode::Mul(args) => {
                for arg in args {
                    let sub = self.transcendental_profile(*arg, depth);
                    profile.merge(&sub);
                }
            }

            ExprNode::Pow { base, exp } => {
                let base_prof = self.transcendental_profile(*base, depth);
                let exp_prof = self.transcendental_profile(*exp, depth);
                profile.merge(&base_prof);
                profile.merge(&exp_prof);

                // Check for fractional powers (algebraic)
                // Note: This is handled by has_radical separately
            }

            ExprNode::Neg(arg) => {
                profile = self.transcendental_profile(*arg, depth);
            }

            ExprNode::Div { num, den } => {
                let num_prof = self.transcendental_profile(*num, depth);
                let den_prof = self.transcendental_profile(*den, depth);
                profile.merge(&num_prof);
                profile.merge(&den_prof);
            }

            ExprNode::Function { id, args } => {
                // Classify the function
                match *id {
                    functions::EXP => {
                        profile.has_exp = true;
                        profile.functions.push(*id);
                    }
                    functions::LN | functions::LOG10 => {
                        profile.has_log = true;
                        profile.functions.push(*id);
                    }
                    functions::SIN | functions::COS | functions::TAN => {
                        profile.has_trig = true;
                        profile.functions.push(*id);
                    }
                    functions::SQRT => {
                        // SQRT is algebraic, not transcendental
                    }
                    _ => {
                        // Unknown function
                        profile.functions.push(*id);
                    }
                }

                // Recurse into arguments
                for arg in args {
                    let sub = self.transcendental_profile(*arg, depth + 1);
                    profile.merge(&sub);
                }
            }
        }

        profile
    }

    /// Checks if expression contains radicals (fractional powers or sqrt).
    fn has_radical(&self, expr: ExprHandle) -> bool {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Integer(_) | ExprNode::Rational(_, _) | ExprNode::Symbol(_) => false,

            ExprNode::Add(args) | ExprNode::Mul(args) => {
                args.iter().any(|h| self.has_radical(*h))
            }

            ExprNode::Pow { base, exp } => {
                // Check if exponent is fractional
                if self.is_fractional_exponent(*exp) {
                    return true;
                }
                self.has_radical(*base) || self.has_radical(*exp)
            }

            ExprNode::Neg(arg) => self.has_radical(*arg),

            ExprNode::Div { num, den } => self.has_radical(*num) || self.has_radical(*den),

            ExprNode::Function { id, args } => {
                if *id == functions::SQRT {
                    return true;
                }
                args.iter().any(|h| self.has_radical(*h))
            }
        }
    }

    /// Checks if an expression represents a fractional exponent.
    fn is_fractional_exponent(&self, expr: ExprHandle) -> bool {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Rational(_, den) => *den != 1,
            ExprNode::Div { .. } => {
                // Any division in exponent likely fractional
                true
            }
            _ => false,
        }
    }

    /// Computes degrees of numerator and denominator for rational expressions.
    /// Returns None if not a rational function.
    fn rational_degrees(&self, expr: ExprHandle) -> Option<(usize, usize)> {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Symbol(id) if *id == self.variable => Some((1, 0)),

            ExprNode::Integer(_) | ExprNode::Rational(_, _) => Some((0, 0)),

            ExprNode::Symbol(_) => Some((0, 0)), // Different variable, treat as constant

            ExprNode::Add(args) => {
                // Degree of sum is max of degrees
                let mut max_num = 0;
                let mut max_den = 0;
                for arg in args {
                    let (n, d) = self.rational_degrees(*arg)?;
                    max_num = max_num.max(n);
                    max_den = max_den.max(d);
                }
                Some((max_num, max_den))
            }

            ExprNode::Mul(args) => {
                // Degree of product is sum of degrees
                let mut sum_num = 0;
                let mut sum_den = 0;
                for arg in args {
                    let (n, d) = self.rational_degrees(*arg)?;
                    sum_num += n;
                    sum_den += d;
                }
                Some((sum_num, sum_den))
            }

            ExprNode::Pow { base, exp } => {
                let (base_num, base_den) = self.rational_degrees(*base)?;
                // Get exponent as integer
                let exp_val = self.expr_to_integer(*exp)?;
                if exp_val >= 0 {
                    let e = exp_val as usize;
                    Some((base_num * e, base_den * e))
                } else {
                    // Negative exponent flips num/den
                    let e = (-exp_val) as usize;
                    Some((base_den * e, base_num * e))
                }
            }

            ExprNode::Neg(arg) => self.rational_degrees(*arg),

            ExprNode::Div { num, den } => {
                let (n_num, n_den) = self.rational_degrees(*num)?;
                let (d_num, d_den) = self.rational_degrees(*den)?;
                // (a/b) / (c/d) = (a*d) / (b*c)
                Some((n_num + d_den, n_den + d_num))
            }

            ExprNode::Function { .. } => None, // Not a rational function
        }
    }

    /// Tries to convert an expression to an integer.
    fn expr_to_integer(&self, expr: ExprHandle) -> Option<i64> {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Integer(n) => Some(*n),
            ExprNode::Neg(arg) => self.expr_to_integer(*arg).map(|n| -n),
            _ => None,
        }
    }

    /// Estimates the degree of the radicand in an algebraic expression.
    fn radicand_degree(&self, expr: ExprHandle) -> usize {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Function { id, args } if *id == functions::SQRT => {
                // Get degree of the argument
                if let Some(arg) = args.first() {
                    self.polynomial_degree(*arg).unwrap_or(0)
                } else {
                    0
                }
            }

            ExprNode::Pow { base, exp } if self.is_fractional_exponent(*exp) => {
                self.polynomial_degree(*base).unwrap_or(0)
            }

            ExprNode::Add(args) | ExprNode::Mul(args) => {
                // Find the maximum radicand degree in subexpressions
                args.iter().map(|h| self.radicand_degree(*h)).max().unwrap_or(0)
            }

            ExprNode::Neg(arg) => self.radicand_degree(*arg),

            ExprNode::Div { num, den } => {
                self.radicand_degree(*num).max(self.radicand_degree(*den))
            }

            _ => 0,
        }
    }

    /// Computes the polynomial degree of an expression in the variable.
    fn polynomial_degree(&self, expr: ExprHandle) -> Option<usize> {
        let node = self.arena.get(expr);
        match node {
            ExprNode::Symbol(id) if *id == self.variable => Some(1),

            ExprNode::Integer(_) | ExprNode::Rational(_, _) => Some(0),

            ExprNode::Symbol(_) => Some(0),

            ExprNode::Add(args) => {
                let mut max_deg = 0;
                for arg in args {
                    max_deg = max_deg.max(self.polynomial_degree(*arg)?);
                }
                Some(max_deg)
            }

            ExprNode::Mul(args) => {
                let mut sum_deg = 0;
                for arg in args {
                    sum_deg += self.polynomial_degree(*arg)?;
                }
                Some(sum_deg)
            }

            ExprNode::Pow { base, exp } => {
                let base_deg = self.polynomial_degree(*base)?;
                let exp_val = self.expr_to_integer(*exp)?;
                if exp_val >= 0 {
                    Some(base_deg * exp_val as usize)
                } else {
                    None // Negative power is not a polynomial
                }
            }

            ExprNode::Neg(arg) => self.polynomial_degree(*arg),

            _ => None,
        }
    }
}

impl TranscendentalProfile {
    /// Merges another profile into this one.
    fn merge(&mut self, other: &TranscendentalProfile) {
        self.has_exp = self.has_exp || other.has_exp;
        self.has_log = self.has_log || other.has_log;
        self.has_trig = self.has_trig || other.has_trig;
        self.has_inverse_trig = self.has_inverse_trig || other.has_inverse_trig;
        self.nesting_depth = self.nesting_depth.max(other.nesting_depth);
        for &f in &other.functions {
            if !self.functions.contains(&f) {
                self.functions.push(f);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;

    #[test]
    fn test_classify_constant() {
        let mut arena = ExprArena::new();
        let x_id = arena.intern_symbol("x");
        let five = arena.integer(5);

        let analyzer = IntegrandAnalyzer::new(&arena, x_id);
        assert_eq!(analyzer.classify(five), IntegrandClass::Constant);
    }

    #[test]
    fn test_classify_polynomial() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");
        let two = arena.integer(2);
        let x_squared = arena.pow(x, two);

        let analyzer = IntegrandAnalyzer::new(&arena, x_id);
        assert_eq!(
            analyzer.classify(x_squared),
            IntegrandClass::Polynomial { degree: 2 }
        );
    }

    #[test]
    fn test_classify_rational() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");
        let one = arena.integer(1);

        // 1/x
        let one_over_x = arena.intern(ExprNode::Div { num: one, den: x });

        let analyzer = IntegrandAnalyzer::new(&arena, x_id);
        let class = analyzer.classify(one_over_x);
        assert!(matches!(class, IntegrandClass::Rational { .. }));
    }

    #[test]
    fn test_classify_transcendental() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");

        // exp(x)
        let exp_x = arena.intern(ExprNode::Function {
            id: functions::EXP,
            args: smallvec![x],
        });

        let analyzer = IntegrandAnalyzer::new(&arena, x_id);
        let class = analyzer.classify(exp_x);
        assert!(matches!(class, IntegrandClass::Transcendental { .. }));
    }

    #[test]
    fn test_classify_algebraic() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let x_id = arena.intern_symbol("x");

        // sqrt(x)
        let sqrt_x = arena.intern(ExprNode::Function {
            id: functions::SQRT,
            args: smallvec![x],
        });

        let analyzer = IntegrandAnalyzer::new(&arena, x_id);
        let class = analyzer.classify(sqrt_x);
        assert!(matches!(class, IntegrandClass::Algebraic { .. }));
    }

    #[test]
    fn test_depends_on_variable() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let y = arena.symbol("y");
        let x_id = arena.intern_symbol("x");

        let analyzer = IntegrandAnalyzer::new(&arena, x_id);

        assert!(analyzer.depends_on_variable(x));
        assert!(!analyzer.depends_on_variable(y));
    }

    #[test]
    fn test_algebraic_genus() {
        assert_eq!(AlgebraicGenus::from_degree(1), AlgebraicGenus::Rational);
        assert_eq!(AlgebraicGenus::from_degree(2), AlgebraicGenus::Rational);
        assert_eq!(AlgebraicGenus::from_degree(3), AlgebraicGenus::Elliptic);
        assert_eq!(AlgebraicGenus::from_degree(4), AlgebraicGenus::Elliptic);
        assert_eq!(AlgebraicGenus::from_degree(5), AlgebraicGenus::Hyperelliptic(2));
    }
}
