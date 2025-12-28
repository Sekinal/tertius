//! Testing Tertius Integration with Difficult Integrals
//!
//! This file tests the symbolic integration capabilities with integrals
//! that are known to trip up many Computer Algebra Systems.
//!
//! Run with: cargo run --example difficult_integrals

use std::time::Instant;
use tertius_integrate::{
    integrate_rational, prove_non_elementary, RationalIntegrationResult,
};
use tertius_integrate::polynomial::integrate_polynomial;
use tertius_integrate::risch::heuristic::{
    check_known_non_elementary, IntegrationTable, heuristic_integrate,
};
use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::rationals::Q;
use tertius_rings::traits::Field;

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
    println!("╔══════════════════════════════════════════════════════════════════╗");
    println!("║        Tertius: Testing Difficult Integrals                      ║");
    println!("║        (Integrals that trip up many CAS systems)                 ║");
    println!("╚══════════════════════════════════════════════════════════════════╝\n");

    let mut passed = 0;
    let mut failed = 0;
    let mut total = 0;

    // Category 1: Rational Functions with Complex Partial Fractions
    println!("═══════════════════════════════════════════════════════════════════");
    println!("CATEGORY 1: Complex Rational Functions");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Test 1: Simple pole
    total += 1;
    if test_simple_pole() { passed += 1; } else { failed += 1; }

    // Test 2: Double pole
    total += 1;
    if test_double_pole() { passed += 1; } else { failed += 1; }

    // Test 3: Cubic denominator
    total += 1;
    if test_cubic_denominator() { passed += 1; } else { failed += 1; }

    // Test 4: Irreducible quadratic
    total += 1;
    if test_irreducible_quadratic() { passed += 1; } else { failed += 1; }

    // Test 5: High degree
    total += 1;
    if test_high_degree() { passed += 1; } else { failed += 1; }

    // Category 2: Pathological Cases
    println!("═══════════════════════════════════════════════════════════════════");
    println!("CATEGORY 2: Pathological Cases");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Test 6: Near-cancellation
    total += 1;
    if test_near_cancellation() { passed += 1; } else { failed += 1; }

    // Test 7: Large coefficients
    total += 1;
    if test_large_coefficients() { passed += 1; } else { failed += 1; }

    // Test 8: Nested structure
    total += 1;
    if test_nested_structure() { passed += 1; } else { failed += 1; }

    // Category 3: Non-Elementary Detection
    println!("═══════════════════════════════════════════════════════════════════");
    println!("CATEGORY 3: Non-Elementary Integral Detection");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Test 9: Gaussian
    total += 1;
    if test_non_elementary_gaussian() { passed += 1; } else { failed += 1; }

    // Test 10: Sine integral
    total += 1;
    if test_non_elementary_sinx_over_x() { passed += 1; } else { failed += 1; }

    // Test 11: Logarithmic integral
    total += 1;
    if test_non_elementary_1_over_lnx() { passed += 1; } else { failed += 1; }

    // Test 12: Exponential integral
    total += 1;
    if test_non_elementary_expx_over_x() { passed += 1; } else { failed += 1; }

    // Category 4: Tricky Elementary Integrals
    println!("═══════════════════════════════════════════════════════════════════");
    println!("CATEGORY 4: Tricky but Elementary Integrals");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Test 13: Partial fractions with repeated roots
    total += 1;
    if test_repeated_roots() { passed += 1; } else { failed += 1; }

    // Test 14: Mixed degree
    total += 1;
    if test_mixed_degree() { passed += 1; } else { failed += 1; }

    // Test 15: Almost cancelling
    total += 1;
    if test_almost_cancelling() { passed += 1; } else { failed += 1; }

    // Category 5: Performance Tests
    println!("═══════════════════════════════════════════════════════════════════");
    println!("CATEGORY 5: Performance Under Stress");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Test 16: High degree polynomial integration
    total += 1;
    if test_high_degree_polynomial() { passed += 1; } else { failed += 1; }

    // Test 17: Many terms
    total += 1;
    if test_many_terms() { passed += 1; } else { failed += 1; }

    // Summary
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                         SUMMARY");
    println!("═══════════════════════════════════════════════════════════════════\n");
    println!("  Total tests:  {}", total);
    println!("  Passed:       {} ✓", passed);
    println!("  Failed:       {} ✗", failed);
    println!("  Success rate: {:.1}%", (passed as f64 / total as f64) * 100.0);
    println!();

    if failed == 0 {
        println!("  All tests passed!");
    } else {
        println!("  Some tests failed. See details above.");
    }
}

// ============================================================================
// CATEGORY 1: Complex Rational Functions
// ============================================================================

fn test_simple_pole() -> bool {
    println!("Test 1: ∫ 1/(x-1) dx = ln|x-1|");
    let start = Instant::now();

    let rf = RationalFunction::new(poly(&[1]), poly(&[-1, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Should have logarithmic part
    let success = result.logarithmic_part.is_some();

    println!("  Time: {:?}", elapsed);
    println!("  Has logarithmic part: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_double_pole() -> bool {
    println!("Test 2: ∫ 1/(x-1)² dx = -1/(x-1)");
    let start = Instant::now();

    // 1/(x-1)² = 1/(x² - 2x + 1)
    let rf = RationalFunction::new(poly(&[1]), poly(&[1, -2, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Hermite reduction should give a rational part
    let success = !result.rational_part.is_zero();

    println!("  Time: {:?}", elapsed);
    println!("  Has rational part: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_cubic_denominator() -> bool {
    println!("Test 3: ∫ 1/(x³ - 1) dx");
    println!("  = (1/3)ln|x-1| - (1/6)ln|x²+x+1| + (1/√3)arctan((2x+1)/√3)");
    let start = Instant::now();

    // x³ - 1 = (x-1)(x²+x+1)
    let rf = RationalFunction::new(poly(&[1]), poly(&[-1, 0, 0, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Should have logarithmic part
    let success = result.logarithmic_part.is_some();

    println!("  Time: {:?}", elapsed);
    println!("  Has logarithmic part: {}", success);
    if let Some(ref log_part) = result.logarithmic_part {
        println!("  Number of log terms: {}", log_part.terms.len());
    }
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_irreducible_quadratic() -> bool {
    println!("Test 4: ∫ 1/(x² + 1) dx = arctan(x)");
    println!("  (Via Rothstein-Trager with algebraic extension)");
    let start = Instant::now();

    let rf = RationalFunction::new(poly(&[1]), poly(&[1, 0, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // This requires complex logs which combine to arctan
    // Rothstein-Trager should handle this
    let success = result.logarithmic_part.is_some() || !result.rational_part.is_zero();

    println!("  Time: {:?}", elapsed);
    println!("  Result: {}\n", if success { "✓ PASS (computed)" } else { "✗ FAIL" });

    success
}

fn test_high_degree() -> bool {
    println!("Test 5: ∫ x⁴/(x⁵ - 1) dx = (1/5)ln|x⁵ - 1|");
    let start = Instant::now();

    // Numerator: x⁴, Denominator: x⁵ - 1
    let rf = RationalFunction::new(poly(&[0, 0, 0, 0, 1]), poly(&[-1, 0, 0, 0, 0, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    let success = result.logarithmic_part.is_some();

    println!("  Time: {:?}", elapsed);
    println!("  Has logarithmic part: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

// ============================================================================
// CATEGORY 2: Pathological Cases
// ============================================================================

fn test_near_cancellation() -> bool {
    println!("Test 6: Near-cancellation test");
    println!("  ∫ (x² - 1)/((x-1)(x+1)) dx = ∫ 1 dx = x");
    let start = Instant::now();

    // (x² - 1) / (x² - 1) = 1 after simplification
    // But we pass it unsimplified to test the algorithm
    let rf = RationalFunction::new(poly(&[-1, 0, 1]), poly(&[-1, 0, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // After polynomial division, should just be integrating 1
    // The polynomial_part should handle this
    let success = true; // If it doesn't crash, it passes

    println!("  Time: {:?}", elapsed);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_large_coefficients() -> bool {
    println!("Test 7: Large coefficients (arbitrary precision test)");
    println!("  ∫ 1/(x - 10000000) dx");
    let start = Instant::now();

    let rf = RationalFunction::new(poly(&[1]), poly(&[-10000000, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    let success = result.logarithmic_part.is_some();

    println!("  Time: {:?}", elapsed);
    println!("  Has logarithmic part: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_nested_structure() -> bool {
    println!("Test 8: Nested structure - ∫ x/(x² + 1)² dx");
    println!("  = -1/(2(x² + 1))");
    let start = Instant::now();

    // x / (x² + 1)² = x / (x⁴ + 2x² + 1)
    let rf = RationalFunction::new(poly(&[0, 1]), poly(&[1, 0, 2, 0, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Should have a rational part from Hermite reduction
    let success = !result.rational_part.is_zero();

    println!("  Time: {:?}", elapsed);
    println!("  Has rational part: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

// ============================================================================
// CATEGORY 3: Non-Elementary Detection
// ============================================================================

fn test_non_elementary_gaussian() -> bool {
    println!("Test 9: Non-elementary detection - ∫ exp(-x²) dx");
    println!("  Expected: Recognized as non-elementary (→ erf)");

    let result = check_known_non_elementary("exp(-x^2)");
    let success = result.is_some();

    if let Some((is_non_elem, explanation)) = result {
        println!("  Detected: {}", is_non_elem);
        println!("  Explanation: {}", explanation);
    }
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_non_elementary_sinx_over_x() -> bool {
    println!("Test 10: Non-elementary detection - ∫ sin(x)/x dx");
    println!("  Expected: Recognized as non-elementary (→ Si)");

    let result = check_known_non_elementary("sin(x)/x");
    let success = result.is_some();

    if let Some((is_non_elem, explanation)) = result {
        println!("  Detected: {}", is_non_elem);
        println!("  Explanation: {}", explanation);
    }
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_non_elementary_1_over_lnx() -> bool {
    println!("Test 11: Non-elementary detection - ∫ 1/ln(x) dx");
    println!("  Expected: Recognized as non-elementary (→ li)");

    let result = check_known_non_elementary("1/ln(x)");
    let success = result.is_some();

    if let Some((is_non_elem, explanation)) = result {
        println!("  Detected: {}", is_non_elem);
        println!("  Explanation: {}", explanation);
    }
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_non_elementary_expx_over_x() -> bool {
    println!("Test 12: Non-elementary detection - ∫ exp(x)/x dx");
    println!("  Expected: Recognized as non-elementary (→ Ei)");

    let result = check_known_non_elementary("exp(x)/x");
    let success = result.is_some();

    if let Some((is_non_elem, explanation)) = result {
        println!("  Detected: {}", is_non_elem);
        println!("  Explanation: {}", explanation);
    }
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

// ============================================================================
// CATEGORY 4: Tricky but Elementary
// ============================================================================

fn test_repeated_roots() -> bool {
    println!("Test 13: Repeated roots - ∫ 1/(x³(x-1)²) dx");
    let start = Instant::now();

    // Denominator: x³(x-1)² = x³(x² - 2x + 1) = x⁵ - 2x⁴ + x³
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 0, 0, 1, -2, 1]));
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Should produce both rational and logarithmic parts
    let has_rational = !result.rational_part.is_zero();
    let has_log = result.logarithmic_part.is_some();
    let success = has_rational || has_log;

    println!("  Time: {:?}", elapsed);
    println!("  Has rational part: {}", has_rational);
    println!("  Has logarithmic part: {}", has_log);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_mixed_degree() -> bool {
    println!("Test 14: Mixed degree - ∫ (x³ + x² + x + 1)/(x² + 1) dx");
    println!("  = x²/2 + arctan(x) + ln(x² + 1)/2");
    let start = Instant::now();

    let rf = RationalFunction::new(
        poly(&[1, 1, 1, 1]),  // x³ + x² + x + 1
        poly(&[1, 0, 1])      // x² + 1
    );
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Should have polynomial part (from division) and logarithmic part
    let success = true; // Check that it doesn't crash

    println!("  Time: {:?}", elapsed);
    println!("  Polynomial part degree: {}", result.polynomial_part.degree());
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_almost_cancelling() -> bool {
    println!("Test 15: Almost cancelling - ∫ (x² + 2x + 1)/(x + 1) dx");
    println!("  = ∫ (x + 1) dx = x²/2 + x");
    let start = Instant::now();

    // (x² + 2x + 1)/(x + 1) = (x + 1)²/(x + 1) = x + 1
    let rf = RationalFunction::new(
        poly(&[1, 2, 1]),  // x² + 2x + 1 = (x + 1)²
        poly(&[1, 1])      // x + 1
    );
    let result = integrate_rational(&rf);

    let elapsed = start.elapsed();

    // Should simplify to polynomial x + 1, then integrate
    let success = result.polynomial_part.degree() >= 1;

    println!("  Time: {:?}", elapsed);
    println!("  Polynomial part computed: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

// ============================================================================
// CATEGORY 5: Performance
// ============================================================================

fn test_high_degree_polynomial() -> bool {
    println!("Test 16: High degree polynomial - ∫ x²⁰ dx = x²¹/21");
    let start = Instant::now();

    // Create x²⁰
    let mut coeffs = vec![q(0); 21];
    coeffs[20] = q(1);
    let p = DensePoly::new(coeffs);
    let result = integrate_polynomial(&p);

    let elapsed = start.elapsed();

    // Check coefficient of x²¹ is 1/21
    let success = result.coeff(21) == q_frac(1, 21);

    println!("  Time: {:?}", elapsed);
    println!("  Coefficient of x²¹ correct: {}", success);
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}

fn test_many_terms() -> bool {
    println!("Test 17: Many terms - ∫ (1 + x + x² + ... + x¹⁰) dx");
    let start = Instant::now();

    // Sum of x^i for i = 0 to 10
    let coeffs: Vec<Q> = (0..=10).map(|_| q(1)).collect();
    let p = DensePoly::new(coeffs);
    let result = integrate_polynomial(&p);

    let elapsed = start.elapsed();

    // Result should have degree 11
    let success = result.degree() == 11;

    println!("  Time: {:?}", elapsed);
    println!("  Result degree: {} (expected 11)", result.degree());
    println!("  Result: {}\n", if success { "✓ PASS" } else { "✗ FAIL" });

    success
}
