//! Exponential and logarithmic simplification rules.

use egg::{rewrite, Rewrite};

use crate::language::TertiusLang;

/// Returns exponential and logarithmic rewrite rules.
pub fn rules() -> Vec<Rewrite<TertiusLang, ()>> {
    vec![
        // exp(0) = 1
        rewrite!("exp-zero"; "(exp 0)" => "1"),

        // ln(1) = 0
        rewrite!("ln-one"; "(ln 1)" => "0"),

        // exp(ln(x)) = x
        rewrite!("exp-ln"; "(exp (ln ?x))" => "?x"),

        // ln(exp(x)) = x
        rewrite!("ln-exp"; "(ln (exp ?x))" => "?x"),

        // ln(x*y) = ln(x) + ln(y)
        rewrite!("ln-mul"; "(ln (* ?x ?y))" => "(+ (ln ?x) (ln ?y))"),
        rewrite!("ln-mul-fold"; "(+ (ln ?x) (ln ?y))" => "(ln (* ?x ?y))"),

        // ln(x/y) = ln(x) - ln(y)
        rewrite!("ln-div"; "(ln (/ ?x ?y))" => "(- (ln ?x) (ln ?y))"),

        // ln(x^n) = n*ln(x)
        rewrite!("ln-pow"; "(ln (^ ?x ?n))" => "(* ?n (ln ?x))"),
        rewrite!("ln-pow-fold"; "(* ?n (ln ?x))" => "(ln (^ ?x ?n))"),

        // exp(a + b) = exp(a) * exp(b)
        rewrite!("exp-add"; "(exp (+ ?a ?b))" => "(* (exp ?a) (exp ?b))"),
        rewrite!("exp-add-fold"; "(* (exp ?a) (exp ?b))" => "(exp (+ ?a ?b))"),

        // exp(a - b) = exp(a) / exp(b)
        rewrite!("exp-sub"; "(exp (- ?a ?b))" => "(/ (exp ?a) (exp ?b))"),

        // exp(n * x) = exp(x)^n
        rewrite!("exp-mul"; "(exp (* ?n ?x))" => "(^ (exp ?x) ?n)"),

        // sqrt(x) = x^(1/2)
        rewrite!("sqrt-pow"; "(sqrt ?x)" => "(^ ?x (/ 1 2))"),
        rewrite!("sqrt-pow-fold"; "(^ ?x (/ 1 2))" => "(sqrt ?x)"),

        // sqrt(x^2) = |x| (needs special handling for sign)
        // For simplification, we often assume positive: sqrt(x^2) = x
        rewrite!("sqrt-sq"; "(sqrt (^ ?x 2))" => "(abs ?x)"),

        // (sqrt(x))^2 = x
        rewrite!("sq-sqrt"; "(^ (sqrt ?x) 2)" => "?x"),

        // sqrt(x*y) = sqrt(x)*sqrt(y)
        rewrite!("sqrt-mul"; "(sqrt (* ?x ?y))" => "(* (sqrt ?x) (sqrt ?y))"),
        rewrite!("sqrt-mul-fold"; "(* (sqrt ?x) (sqrt ?y))" => "(sqrt (* ?x ?y))"),

        // ln(sqrt(x)) = ln(x)/2
        rewrite!("ln-sqrt"; "(ln (sqrt ?x))" => "(/ (ln ?x) 2)"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use egg::Runner;

    #[test]
    fn test_exp_ln() {
        let rules = rules();
        let start = "(exp (ln x))".parse().unwrap();
        let runner = Runner::default().with_expr(&start).run(&rules);
        let extractor = egg::Extractor::new(&runner.egraph, egg::AstSize);
        let (_, best) = extractor.find_best(runner.roots[0]);
        assert_eq!(best.to_string(), "x");
    }

    #[test]
    fn test_ln_exp() {
        let rules = rules();
        let start = "(ln (exp x))".parse().unwrap();
        let runner = Runner::default().with_expr(&start).run(&rules);
        let extractor = egg::Extractor::new(&runner.egraph, egg::AstSize);
        let (_, best) = extractor.find_best(runner.roots[0]);
        assert_eq!(best.to_string(), "x");
    }
}
