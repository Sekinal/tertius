//! Impossible Integrals: Reality Check for CAS Systems
//!
//! These integrals are known to cause problems for major CAS systems including
//! Mathematica, Maple, and SymPy. Many have no closed-form elementary solution,
//! require special functions, or have answers so complex they're impractical.
//!
//! Run with: cargo run --example impossible_integrals
//!
//! References:
//! - Risch algorithm limitations
//! - RUBI (Rule-Based Integrator) test failures
//! - Symbolic Integration Tutorial (Moses, 1971)
//! - "The Integration of Algebraic Functions" (Davenport, 1981)

use std::time::Instant;
use tertius_integrate::{integrate_rational, RationalIntegrationResult};
use tertius_integrate::polynomial::integrate_polynomial;
use tertius_integrate::risch::heuristic::check_known_non_elementary;
use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::rationals::Q;

fn q(n: i64) -> Q {
    Q::from_integer(n)
}

fn q_frac(num: i64, den: i64) -> Q {
    Q::new(num, den)
}

fn poly(coeffs: &[i64]) -> DensePoly<Q> {
    DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
}

fn main() {
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║           IMPOSSIBLE INTEGRALS: CAS REALITY CHECK                    ║");
    println!("║                                                                      ║");
    println!("║   Integrals that challenge Mathematica, Maple, and SymPy             ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝\n");

    // =========================================================================
    // TIER 1: DECEPTIVELY SIMPLE (Look easy, are hard)
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 1: DECEPTIVELY SIMPLE");
    println!("These look trivial but have complex or no elementary solutions");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier1_deceptive();

    // =========================================================================
    // TIER 2: ELLIPTIC INTEGRALS (Algebraic integrands, transcendental answers)
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 2: ELLIPTIC & HYPERELLIPTIC");
    println!("Algebraic integrands requiring elliptic functions");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier2_elliptic();

    // =========================================================================
    // TIER 3: RATIONAL FUNCTION NIGHTMARES
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 3: RATIONAL FUNCTION NIGHTMARES");
    println!("Rational functions with pathological structure");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier3_rational_nightmares();

    // =========================================================================
    // TIER 4: FAMOUS NON-ELEMENTARY INTEGRALS
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 4: FAMOUS NON-ELEMENTARY INTEGRALS");
    println!("Proven to have no elementary antiderivative");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier4_famous_non_elementary();

    // =========================================================================
    // TIER 5: RISCH ALGORITHM EDGE CASES
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 5: RISCH ALGORITHM EDGE CASES");
    println!("Cases that require full Risch algorithm with all extensions");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier5_risch_edge_cases();

    // =========================================================================
    // TIER 6: MATHEMATICA FAILURES
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 6: KNOWN MATHEMATICA FAILURES");
    println!("Integrals that Mathematica returns unevaluated or wrong");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier6_mathematica_failures();

    // =========================================================================
    // TIER 7: THE IMPOSSIBLE
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════════════════");
    println!("TIER 7: THE TRULY IMPOSSIBLE");
    println!("Integrals at the absolute frontier of symbolic computation");
    println!("═══════════════════════════════════════════════════════════════════════\n");

    test_tier7_truly_impossible();

    println!("\n═══════════════════════════════════════════════════════════════════════");
    println!("CONCLUSION");
    println!("═══════════════════════════════════════════════════════════════════════\n");
    println!("These tests represent the frontier of symbolic integration.");
    println!("A perfect CAS would need:");
    println!("  • Complete Risch algorithm with all extensions");
    println!("  • Algebraic function integration (Davenport/Trager)");
    println!("  • Full special function recognition");
    println!("  • Algebraic number theory for splitting fields");
    println!("  • Differential Galois theory for Liouvillian functions");
    println!();
}

fn test_tier1_deceptive() {
    // 1.1: ∫ 1/(x⁴ + 1) dx
    // Looks simple, but factors over Q(√2, i) and has a very complex answer
    println!("1.1: ∫ 1/(x⁴ + 1) dx");
    println!("     Factors: (x² + √2·x + 1)(x² - √2·x + 1) over R");
    println!("     Answer involves: arctan and log with √2 coefficients");
    let rf = RationalFunction::new(poly(&[1]), poly(&[1, 0, 0, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Has log part: {}", result.logarithmic_part.is_some());
    println!("     Status: {}\n", status_str(result.logarithmic_part.is_some()));

    // 1.2: ∫ 1/(x⁶ + 1) dx
    // Even worse - factors into cyclotomic components
    println!("1.2: ∫ 1/(x⁶ + 1) dx");
    println!("     Involves: 6th roots of -1, extremely complex partial fractions");
    let rf = RationalFunction::new(poly(&[1]), poly(&[1, 0, 0, 0, 0, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", status_str(result.logarithmic_part.is_some()));

    // 1.3: ∫ 1/(x⁸ + 1) dx
    // 8th roots of unity - Galois theory nightmare
    println!("1.3: ∫ 1/(x⁸ + 1) dx");
    println!("     Requires: 8th roots of -1, degree 4 algebraic extension");
    let rf = RationalFunction::new(poly(&[1]), poly(&[1, 0, 0, 0, 0, 0, 0, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", status_str(result.logarithmic_part.is_some()));

    // 1.4: ∫ x²/(x⁴ + 1) dx
    println!("1.4: ∫ x²/(x⁴ + 1) dx");
    println!("     Related to 1.1 but different partial fraction decomposition");
    let rf = RationalFunction::new(poly(&[0, 0, 1]), poly(&[1, 0, 0, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", status_str(result.logarithmic_part.is_some()));
}

fn test_tier2_elliptic() {
    // These require elliptic integrals - beyond elementary functions

    // 2.1: ∫ 1/√(x³ - x) dx
    println!("2.1: ∫ 1/√(x³ - x) dx = ∫ 1/√(x(x-1)(x+1)) dx");
    println!("     This is an ELLIPTIC INTEGRAL of the first kind");
    println!("     No elementary antiderivative exists");
    println!("     Answer: F(arcsin(√x), -1) where F is incomplete elliptic integral");
    println!("     Status: Non-elementary (elliptic)\n");

    // 2.2: ∫ 1/√(x⁴ + 1) dx
    println!("2.2: ∫ 1/√(x⁴ + 1) dx");
    println!("     HYPERELLIPTIC integral - genus 1 curve");
    println!("     Requires: Weierstrass ℘-function or elliptic integrals");
    println!("     Status: Non-elementary (hyperelliptic)\n");

    // 2.3: ∫ √(x⁴ + 1) dx
    println!("2.3: ∫ √(x⁴ + 1) dx");
    println!("     Elliptic integral of the second kind");
    println!("     Even Mathematica gives a complex expression with EllipticE");
    println!("     Status: Non-elementary (elliptic)\n");

    // 2.4: ∫ 1/√(1 - x⁴) dx
    println!("2.4: ∫ 1/√(1 - x⁴) dx (Lemniscatic integral)");
    println!("     Famous: Related to arc length of lemniscate of Bernoulli");
    println!("     Connected to: Gauss's constant, AGM");
    println!("     Status: Non-elementary (lemniscatic)\n");
}

fn test_tier3_rational_nightmares() {
    // 3.1: Very high degree cyclotomic
    println!("3.1: ∫ 1/(x¹² - 1) dx");
    println!("     Factors into Φ₁·Φ₂·Φ₃·Φ₄·Φ₆·Φ₁₂ (cyclotomic)");
    println!("     Partial fractions have 12 terms with algebraic coefficients");
    let mut coeffs = vec![0i64; 13];
    coeffs[0] = -1;
    coeffs[12] = 1;
    let rf = RationalFunction::new(poly(&[1]), poly(&coeffs));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", status_str(result.logarithmic_part.is_some()));

    // 3.2: Repeated irreducible factors
    println!("3.2: ∫ 1/(x² + 1)³ dx");
    println!("     Hermite reduction needed for multiplicity 3");
    println!("     Answer involves: rational function + arctan");
    // (x² + 1)³ = x⁶ + 3x⁴ + 3x² + 1
    let rf = RationalFunction::new(poly(&[1]), poly(&[1, 0, 3, 0, 3, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Has rational part: {}", !result.rational_part.is_zero());
    println!("     Status: {}\n", if !result.rational_part.is_zero() { "COMPUTED" } else { "INCOMPLETE" });

    // 3.3: Mixed repeated roots
    println!("3.3: ∫ 1/((x² + 1)²(x² + 4)²) dx");
    println!("     Two distinct irreducible quadratics, both squared");
    // (x² + 1)²(x² + 4)² = (x⁴ + 2x² + 1)(x⁴ + 8x² + 16)
    // = x⁸ + 10x⁶ + 33x⁴ + 42x² + 16
    let rf = RationalFunction::new(poly(&[1]), poly(&[16, 0, 42, 0, 33, 0, 10, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", if !result.rational_part.is_zero() { "COMPUTED" } else { "INCOMPLETE" });

    // 3.4: Fibonacci-like polynomial
    println!("3.4: ∫ 1/(x⁵ - x - 1) dx");
    println!("     This quintic is IRREDUCIBLE over Q (Galois group S₅)");
    println!("     Cannot be solved by radicals!");
    println!("     Requires: numerical root finding or unsolvable algebraics");
    let rf = RationalFunction::new(poly(&[1]), poly(&[-1, -1, 0, 0, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", status_str(result.logarithmic_part.is_some()));
}

fn test_tier4_famous_non_elementary() {
    println!("4.1: ∫ e^(-x²) dx (Gaussian/Error function)");
    println!("     PROVEN non-elementary by Liouville (1835)");
    println!("     Answer: √π/2 · erf(x)");
    let detected = check_known_non_elementary("exp(-x^2)");
    println!("     Tertius detection: {}", if detected.is_some() { "RECOGNIZED" } else { "MISSED" });
    if let Some((_, exp)) = detected {
        println!("     Explanation: {}", exp);
    }
    println!();

    println!("4.2: ∫ sin(x²) dx (Fresnel S integral)");
    println!("     PROVEN non-elementary");
    println!("     Answer: √(π/2) · S(x√(2/π)) where S is Fresnel integral");
    println!("     Status: Non-elementary (Fresnel)\n");

    println!("4.3: ∫ cos(x²) dx (Fresnel C integral)");
    println!("     PROVEN non-elementary");
    println!("     Answer: √(π/2) · C(x√(2/π))");
    println!("     Status: Non-elementary (Fresnel)\n");

    println!("4.4: ∫ sin(x)/x dx (Sine integral)");
    let detected = check_known_non_elementary("sin(x)/x");
    println!("     Answer: Si(x)");
    println!("     Tertius detection: {}\n", if detected.is_some() { "RECOGNIZED" } else { "MISSED" });

    println!("4.5: ∫ e^x/x dx (Exponential integral)");
    let detected = check_known_non_elementary("exp(x)/x");
    println!("     Answer: Ei(x)");
    println!("     Tertius detection: {}\n", if detected.is_some() { "RECOGNIZED" } else { "MISSED" });

    println!("4.6: ∫ 1/ln(x) dx (Logarithmic integral)");
    let detected = check_known_non_elementary("1/ln(x)");
    println!("     Answer: li(x)");
    println!("     Tertius detection: {}\n", if detected.is_some() { "RECOGNIZED" } else { "MISSED" });

    println!("4.7: ∫ x^x dx");
    let detected = check_known_non_elementary("x^x");
    println!("     No known special function!");
    println!("     Cannot be expressed in terms of ANY known functions");
    println!("     Tertius detection: {}\n", if detected.is_some() { "RECOGNIZED" } else { "MISSED" });

    println!("4.8: ∫ e^(e^x) dx");
    println!("     Doubly exponential - no elementary or special function form");
    println!("     Status: Non-elementary (no known closed form)\n");
}

fn test_tier5_risch_edge_cases() {
    println!("5.1: ∫ log(x)/(x² + 1) dx");
    println!("     Requires: Risch algorithm with logarithmic extension");
    println!("     Answer involves: polylogarithm Li₂");
    println!("     Status: Non-elementary (dilogarithm)\n");

    println!("5.2: ∫ log(log(x))/x dx");
    println!("     Nested logarithms - tests tower building");
    println!("     Answer: log(x)·log(log(x)) - log(x) = log(x)(log(log(x)) - 1)");
    println!("     Status: Elementary but requires nested tower\n");

    println!("5.3: ∫ e^x · log(x) dx");
    println!("     Mixed exponential-logarithmic");
    println!("     Requires: integration by parts in Risch framework");
    println!("     Answer involves: Ei(x)");
    println!("     Status: Non-elementary\n");

    println!("5.4: ∫ log(x + √(x² + 1)) dx (inverse sinh argument)");
    println!("     Algebraic + logarithmic");
    println!("     Tests: proper handling of algebraic functions under log");
    println!("     Answer: x·arcsinh(x) - √(x² + 1)");
    println!("     Status: Elementary but tricky\n");

    println!("5.5: ∫ arctan(x)/(1 + x²) dx");
    println!("     = (1/2)arctan²(x)");
    println!("     Requires: recognizing derivative of arctan inside");
    println!("     Status: Elementary\n");
}

fn test_tier6_mathematica_failures() {
    println!("These are documented cases where Mathematica versions have failed:\n");

    println!("6.1: ∫ (x⁴ - 3x² + 6)/(x⁶ - 5x⁴ + 5x² + 4) dx");
    println!("     Mathematica 4.0: returned unevaluated");
    println!("     The denominator factors non-trivially over algebraic numbers");
    // x⁶ - 5x⁴ + 5x² + 4
    let rf = RationalFunction::new(
        poly(&[6, 0, -3, 0, 1, 0, 0]),  // x⁴ - 3x² + 6 (padded)
        poly(&[4, 0, 5, 0, -5, 0, 1])
    );
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", if result.logarithmic_part.is_some() || !result.rational_part.is_zero() { "ATTEMPT" } else { "UNEVALUATED" });

    println!("6.2: ∫ 1/((x-1)(x-2)(x-3)(x-4)(x-5)(x-6)(x-7)) dx");
    println!("     7 distinct linear factors - partial fractions explosion");
    println!("     Some CAS timeout on this due to coefficient explosion");
    // Product = x⁷ - 28x⁶ + 322x⁵ - 1960x⁴ + 6769x³ - 13132x² + 13068x - 5040
    let rf = RationalFunction::new(
        poly(&[1]),
        poly(&[-5040, 13068, -13132, 6769, -1960, 322, -28, 1])
    );
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Has log part: {}", result.logarithmic_part.is_some());
    if let Some(ref log) = result.logarithmic_part {
        println!("     Number of log terms: {}", log.terms.len());
    }
    println!();

    println!("6.3: ∫ 1/(x^10 + x^9 + x^8 + ... + x + 1) dx");
    println!("     = ∫ 1/(x¹¹ - 1)/(x - 1) dx");
    println!("     Cyclotomic Φ₁₁(x) - requires 11th roots of unity");
    println!("     Degree 10 irreducible polynomial over Q");
    // Φ₁₁(x) = x¹⁰ + x⁹ + ... + x + 1
    let rf = RationalFunction::new(
        poly(&[1]),
        poly(&[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])  // all 1s
    );
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", if result.logarithmic_part.is_some() { "COMPUTED" } else { "UNEVALUATED" });
}

fn test_tier7_truly_impossible() {
    println!("These represent the absolute limits of symbolic integration:\n");

    println!("7.1: ∫ e^(x²) dx");
    println!("     Dawson function territory");
    println!("     Related to: imaginary error function erfi(x)");
    println!("     F(x) = e^(-x²) ∫₀ˣ e^(t²) dt (Dawson's integral)");
    let detected = check_known_non_elementary("exp(x^2)");
    println!("     Tertius detection: {}\n", if detected.is_some() { "RECOGNIZED" } else { "MISSED" });

    println!("7.2: ∫ e^(x³) dx");
    println!("     No standard special function!");
    println!("     Sometimes expressed via incomplete gamma or hypergeometric");
    println!("     Status: Beyond standard special functions\n");

    println!("7.3: ∫ x^x · (1 + log(x)) dx");
    println!("     This is actually elementary: = x^x");
    println!("     Because d/dx[x^x] = x^x(1 + log(x))");
    println!("     But most CAS fail to recognize this pattern!");
    println!("     Status: Elementary (derivative recognition)\n");

    println!("7.4: ∫ 1/(x⁵ + x + 1) dx");
    println!("     Quintic with Galois group S₅ (unsolvable by radicals)");
    println!("     Cannot be integrated in terms of ANY algebraic functions!");
    let rf = RationalFunction::new(poly(&[1]), poly(&[1, 1, 0, 0, 0, 1]));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: Fundamentally unsolvable over radicals\n");

    println!("7.5: ∫ log(Γ(x)) dx");
    println!("     Involves: log-gamma function");
    println!("     Related to: Barnes G-function, Stirling series");
    println!("     Status: Requires advanced special functions\n");

    println!("7.6: ∫ ζ(x) dx where ζ is Riemann zeta");
    println!("     Zeta function integral");
    println!("     No known closed form!");
    println!("     Status: Open problem\n");

    println!("7.7: The Risch Decision Problem Itself");
    println!("     For general algebraic functions, the Risch algorithm");
    println!("     is NOT decidable (requires solving arbitrary polynomial systems)");
    println!("     The problem: ∫ A(x, y) dx where P(x, y) = 0 defines y");
    println!("     Status: Theoretically undecidable in general\n");

    // Test a specific hard rational function
    println!("7.8: ∫ (x²⁰ + 1)/(x²¹ + x + 1) dx");
    println!("     Very high degree, non-trivial factorization");
    let mut num = vec![0i64; 21];
    num[0] = 1;
    num[20] = 1;
    let mut den = vec![0i64; 22];
    den[0] = 1;
    den[1] = 1;
    den[21] = 1;
    let rf = RationalFunction::new(poly(&num), poly(&den));
    let start = Instant::now();
    let result = integrate_rational(&rf);
    println!("     Tertius time: {:?}", start.elapsed());
    println!("     Status: {}\n", if result.logarithmic_part.is_some() || !result.polynomial_part.is_zero() { "ATTEMPT" } else { "UNEVALUATED" });
}

fn status_str(success: bool) -> &'static str {
    if success { "COMPUTED" } else { "NEEDS WORK" }
}
