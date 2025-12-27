//! Basic algebraic simplification rules.

use egg::{rewrite, Rewrite};

use crate::language::TertiusLang;

/// Returns basic arithmetic rewrite rules.
pub fn rules() -> Vec<Rewrite<TertiusLang, ()>> {
    vec![
        // Additive identity
        rewrite!("add-zero-l"; "(+ 0 ?a)" => "?a"),
        rewrite!("add-zero-r"; "(+ ?a 0)" => "?a"),

        // Multiplicative identity
        rewrite!("mul-one-l"; "(* 1 ?a)" => "?a"),
        rewrite!("mul-one-r"; "(* ?a 1)" => "?a"),

        // Multiplicative zero
        rewrite!("mul-zero-l"; "(* 0 ?a)" => "0"),
        rewrite!("mul-zero-r"; "(* ?a 0)" => "0"),

        // Double negation
        rewrite!("neg-neg"; "(neg (neg ?a))" => "?a"),

        // Subtraction as addition of negation
        rewrite!("sub-to-add"; "(- ?a ?b)" => "(+ ?a (neg ?b))"),

        // Commutativity
        rewrite!("add-comm"; "(+ ?a ?b)" => "(+ ?b ?a)"),
        rewrite!("mul-comm"; "(* ?a ?b)" => "(* ?b ?a)"),

        // Associativity
        rewrite!("add-assoc-l"; "(+ (+ ?a ?b) ?c)" => "(+ ?a (+ ?b ?c))"),
        rewrite!("add-assoc-r"; "(+ ?a (+ ?b ?c))" => "(+ (+ ?a ?b) ?c)"),
        rewrite!("mul-assoc-l"; "(* (* ?a ?b) ?c)" => "(* ?a (* ?b ?c))"),
        rewrite!("mul-assoc-r"; "(* ?a (* ?b ?c))" => "(* (* ?a ?b) ?c)"),

        // Distributivity
        rewrite!("dist-l"; "(* ?a (+ ?b ?c))" => "(+ (* ?a ?b) (* ?a ?c))"),
        rewrite!("dist-r"; "(* (+ ?a ?b) ?c)" => "(+ (* ?a ?c) (* ?b ?c))"),

        // Factor (reverse of distributivity)
        rewrite!("factor-l"; "(+ (* ?a ?b) (* ?a ?c))" => "(* ?a (+ ?b ?c))"),

        // Power rules
        rewrite!("pow-zero"; "(^ ?a 0)" => "1"),
        rewrite!("pow-one"; "(^ ?a 1)" => "?a"),
        rewrite!("pow-neg-one"; "(^ ?a (neg 1))" => "(/ 1 ?a)"),

        // Power of product
        rewrite!("pow-mul"; "(^ (* ?a ?b) ?n)" => "(* (^ ?a ?n) (^ ?b ?n))"),

        // Power of power
        rewrite!("pow-pow"; "(^ (^ ?a ?m) ?n)" => "(^ ?a (* ?m ?n))"),

        // Product of powers (same base)
        rewrite!("mul-pow"; "(* (^ ?a ?m) (^ ?a ?n))" => "(^ ?a (+ ?m ?n))"),

        // Division rules
        rewrite!("div-one"; "(/ ?a 1)" => "?a"),
        rewrite!("div-self"; "(/ ?a ?a)" => "1"),
        rewrite!("div-neg"; "(/ (neg ?a) ?b)" => "(neg (/ ?a ?b))"),

        // a + a = 2a
        rewrite!("add-same"; "(+ ?a ?a)" => "(* 2 ?a)"),

        // a * a = a^2
        rewrite!("mul-same"; "(* ?a ?a)" => "(^ ?a 2)"),

        // Negative multiplication
        rewrite!("neg-mul-l"; "(* (neg ?a) ?b)" => "(neg (* ?a ?b))"),
        rewrite!("neg-mul-r"; "(* ?a (neg ?b))" => "(neg (* ?a ?b))"),
        rewrite!("neg-mul-both"; "(* (neg ?a) (neg ?b))" => "(* ?a ?b)"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use egg::{EGraph, Runner};

    #[test]
    fn test_add_zero() {
        let rules = rules();
        let start = "(+ x 0)".parse().unwrap();
        let runner = Runner::default().with_expr(&start).run(&rules);
        let extractor = egg::Extractor::new(&runner.egraph, egg::AstSize);
        let (_, best) = extractor.find_best(runner.roots[0]);
        assert_eq!(best.to_string(), "x");
    }

    #[test]
    fn test_mul_zero() {
        let rules = rules();
        let start = "(* x 0)".parse().unwrap();
        let runner = Runner::default().with_expr(&start).run(&rules);
        let extractor = egg::Extractor::new(&runner.egraph, egg::AstSize);
        let (_, best) = extractor.find_best(runner.roots[0]);
        assert_eq!(best.to_string(), "0");
    }
}
