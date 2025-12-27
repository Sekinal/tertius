//! Trigonometric simplification rules.

use egg::{rewrite, Rewrite};

use crate::language::TertiusLang;

/// Returns trigonometric rewrite rules.
pub fn rules() -> Vec<Rewrite<TertiusLang, ()>> {
    vec![
        // Pythagorean identity: sin²(x) + cos²(x) = 1
        rewrite!("pythag"; "(+ (^ (sin ?x) 2) (^ (cos ?x) 2))" => "1"),
        rewrite!("pythag-r"; "(+ (^ (cos ?x) 2) (^ (sin ?x) 2))" => "1"),

        // tan = sin/cos
        rewrite!("tan-def"; "(tan ?x)" => "(/ (sin ?x) (cos ?x))"),
        rewrite!("tan-fold"; "(/ (sin ?x) (cos ?x))" => "(tan ?x)"),

        // sin(0) = 0, cos(0) = 1
        rewrite!("sin-zero"; "(sin 0)" => "0"),
        rewrite!("cos-zero"; "(cos 0)" => "1"),

        // sin(-x) = -sin(x) (odd function)
        rewrite!("sin-neg"; "(sin (neg ?x))" => "(neg (sin ?x))"),

        // cos(-x) = cos(x) (even function)
        rewrite!("cos-neg"; "(cos (neg ?x))" => "(cos ?x)"),

        // tan(-x) = -tan(x) (odd function)
        rewrite!("tan-neg"; "(tan (neg ?x))" => "(neg (tan ?x))"),

        // Double angle formulas
        // sin(2x) = 2*sin(x)*cos(x)
        rewrite!("sin-double"; "(sin (* 2 ?x))" => "(* 2 (* (sin ?x) (cos ?x)))"),

        // cos(2x) = cos²(x) - sin²(x)
        rewrite!("cos-double"; "(cos (* 2 ?x))" => "(- (^ (cos ?x) 2) (^ (sin ?x) 2))"),

        // Sum-to-product (for special cases)
        // sin²(x) = (1 - cos(2x))/2
        rewrite!("sin-sq"; "(^ (sin ?x) 2)" => "(/ (- 1 (cos (* 2 ?x))) 2)"),

        // cos²(x) = (1 + cos(2x))/2
        rewrite!("cos-sq"; "(^ (cos ?x) 2)" => "(/ (+ 1 (cos (* 2 ?x))) 2)"),

        // Inverse function compositions
        rewrite!("sin-asin"; "(sin (asin ?x))" => "?x"),
        rewrite!("cos-acos"; "(cos (acos ?x))" => "?x"),
        rewrite!("tan-atan"; "(tan (atan ?x))" => "?x"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use egg::Runner;

    #[test]
    fn test_pythag() {
        let rules = rules();
        let start = "(+ (^ (sin x) 2) (^ (cos x) 2))".parse().unwrap();
        let runner = Runner::default().with_expr(&start).run(&rules);
        let extractor = egg::Extractor::new(&runner.egraph, egg::AstSize);
        let (_, best) = extractor.find_best(runner.roots[0]);
        assert_eq!(best.to_string(), "1");
    }
}
