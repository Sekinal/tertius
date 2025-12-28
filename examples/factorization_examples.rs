//! Multivariate Polynomial Factorization Examples
//!
//! Demonstrates Tertius's capabilities for factoring polynomials in multiple variables
//! using Lecerf's algorithm and Van Hoeij for univariate cases.
//!
//! Run with: cargo run --example factorization_examples

use tertius_factor::{
    factor_multivariate, lecerf_factor_multivariate, van_hoeij_factor, MultivariateFactorResult,
};
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_poly::monomial::PackedMonomial;
use tertius_poly::ordering::MonomialOrder;
use tertius_poly::sparse::SparsePoly;
use tertius_rings::integers::Z;
use std::time::Instant;

// Helper to create integers
fn z(n: i64) -> Z {
    Z(Integer::new(n))
}

// Helper to create dense polynomials over Z
fn poly_z(coeffs: &[i64]) -> DensePoly<Z> {
    DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
}

// Helper to create monomials
fn mono(exps: &[u32]) -> PackedMonomial {
    PackedMonomial::from_exponents(exps)
}

// Helper to create sparse polynomials
fn sparse_poly(terms: &[(i64, &[u32])], num_vars: usize) -> SparsePoly<Z> {
    let t: Vec<_> = terms.iter().map(|(c, e)| (mono(e), z(*c))).collect();
    SparsePoly::new(t, num_vars, MonomialOrder::Grevlex)
}

fn main() {
    println!("╔════════════════════════════════════════════════════════════════════╗");
    println!("║     Tertius: Multivariate Polynomial Factorization Examples        ║");
    println!("╚════════════════════════════════════════════════════════════════════╝\n");

    example_1_univariate_factorization();
    example_2_bivariate_factorization();
    example_3_trivariate_factorization();
    example_4_dispatch_function();
    example_5_special_polynomials();
    example_6_performance_comparison();
}

/// Example 1: Univariate Factorization (Van Hoeij Algorithm)
fn example_1_univariate_factorization() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 1: Univariate Factorization (Van Hoeij + LLL)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Factor x⁶ - 1 = (x-1)(x+1)(x²+x+1)(x²-x+1)
    let p1 = poly_z(&[-1, 0, 0, 0, 0, 0, 1]);
    println!("  Factoring: x⁶ - 1 (cyclotomic factorization)");

    let start = Instant::now();
    let result = van_hoeij_factor(&p1);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} irreducible factors:", result.factors.len());
    for (i, f) in result.factors.iter().enumerate() {
        let coeffs: Vec<i64> = f.coeffs().iter().map(|c| {
            // Convert to i64 for display
            let s = format!("{:?}", c.0);
            s.parse().unwrap_or(0)
        }).collect();
        println!("    Factor {}: {:?}", i + 1, coeffs);
    }
    println!();

    // Factor x⁸ - 1
    let p2 = poly_z(&[-1, 0, 0, 0, 0, 0, 0, 0, 1]);
    println!("  Factoring: x⁸ - 1");

    let start = Instant::now();
    let result = van_hoeij_factor(&p2);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} irreducible factors", result.factors.len());
    println!();

    // Factor Swinnerton-Dyer polynomial (challenging for many CAS)
    // S_2 = (x - √2 - √3)(x - √2 + √3)(x + √2 - √3)(x + √2 + √3)
    //     = x⁴ - 10x² + 1
    let swinnerton = poly_z(&[1, 0, -10, 0, 1]);
    println!("  Factoring: x⁴ - 10x² + 1 (Swinnerton-Dyer S₂)");
    println!("  This polynomial is irreducible over Z but factors over Q(√2, √3)");

    let start = Instant::now();
    let result = van_hoeij_factor(&swinnerton);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Irreducible: {}", result.factors.len() == 1);
    println!();
}

/// Example 2: Bivariate Factorization
fn example_2_bivariate_factorization() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 2: Bivariate Factorization");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // f = (x + 1)(x + 2) = x² + 3x + 2 in Z[x, y]
    let f1 = sparse_poly(&[(1, &[2, 0]), (3, &[1, 0]), (2, &[0, 0])], 2);
    println!("  Factoring: x² + 3x + 2 in Z[x, y]");
    println!("  Expected: (x + 1)(x + 2)");

    let start = Instant::now();
    let result = factor_multivariate(&f1);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} factors", result.factors.len());
    print_factors(&result);
    println!();

    // f = x² - y² = (x - y)(x + y)
    let f2 = sparse_poly(&[(1, &[2, 0]), (-1, &[0, 2])], 2);
    println!("  Factoring: x² - y² (difference of squares)");
    println!("  Expected: (x - y)(x + y)");

    let start = Instant::now();
    let result = factor_multivariate(&f2);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} factors", result.factors.len());
    print_factors(&result);
    println!();

    // f = x² + 2xy + y² = (x + y)²
    let f3 = sparse_poly(&[(1, &[2, 0]), (2, &[1, 1]), (1, &[0, 2])], 2);
    println!("  Factoring: x² + 2xy + y² (perfect square)");
    println!("  Expected: (x + y)²");

    let start = Instant::now();
    let result = factor_multivariate(&f3);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} factors", result.factors.len());
    print_factors(&result);
    println!();
}

/// Example 3: Trivariate Factorization (Lecerf's Algorithm)
fn example_3_trivariate_factorization() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 3: Trivariate Factorization (Lecerf's Algorithm)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // f = x² + 2x + 1 = (x + 1)² in Z[x, y, z]
    let f1 = sparse_poly(&[(1, &[2, 0, 0]), (2, &[1, 0, 0]), (1, &[0, 0, 0])], 3);
    println!("  Factoring: x² + 2x + 1 in Z[x, y, z]");
    println!("  Expected: (x + 1)²");

    let start = Instant::now();
    let result = lecerf_factor_multivariate(&f1);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} factors", result.factors.len());
    println!("  Lifting success: {}", result.stats.lifting_success);
    println!();

    // f = (x + 1)(x - 1) = x² - 1 in Z[x, y, z]
    let f2 = sparse_poly(&[(1, &[2, 0, 0]), (-1, &[0, 0, 0])], 3);
    println!("  Factoring: x² - 1 in Z[x, y, z]");
    println!("  Expected: (x + 1)(x - 1)");

    let start = Instant::now();
    let result = lecerf_factor_multivariate(&f2);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Found {} factors", result.factors.len());
    println!();

    // Irreducible trivariate: x² + y² + z² (sum of three squares)
    let f3 = sparse_poly(&[(1, &[2, 0, 0]), (1, &[0, 2, 0]), (1, &[0, 0, 2])], 3);
    println!("  Factoring: x² + y² + z² (irreducible over Z)");

    let start = Instant::now();
    let result = lecerf_factor_multivariate(&f3);
    let elapsed = start.elapsed();

    println!("  Time: {:?}", elapsed);
    println!("  Irreducible: {}", result.factors.len() == 1);
    println!();
}

/// Example 4: Unified Dispatch Function
fn example_4_dispatch_function() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 4: Unified factor_multivariate() Dispatch");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    println!("  The factor_multivariate() function automatically dispatches to:");
    println!("    • Van Hoeij for univariate polynomials");
    println!("    • Bivariate Hensel for 2 variables");
    println!("    • Lecerf's algorithm for 3+ variables\n");

    // Univariate (dispatches to Van Hoeij)
    let uni = sparse_poly(&[(1, &[2]), (3, &[1]), (2, &[0])], 1);
    println!("  1 variable: x² + 3x + 2");
    let result = factor_multivariate(&uni);
    println!("    → {} factors (Van Hoeij)", result.factors.len());

    // Bivariate (dispatches to bivariate Hensel)
    let biv = sparse_poly(&[(1, &[2, 0]), (3, &[1, 0]), (2, &[0, 0])], 2);
    println!("  2 variables: x² + 3x + 2 in Z[x,y]");
    let result = factor_multivariate(&biv);
    println!("    → {} factors (Bivariate Hensel)", result.factors.len());

    // Trivariate (dispatches to Lecerf)
    let triv = sparse_poly(&[(1, &[2, 0, 0]), (3, &[1, 0, 0]), (2, &[0, 0, 0])], 3);
    println!("  3 variables: x² + 3x + 2 in Z[x,y,z]");
    let result = factor_multivariate(&triv);
    println!("    → {} factors (Lecerf)", result.factors.len());
    println!();
}

/// Example 5: Special Polynomials
fn example_5_special_polynomials() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 5: Special Polynomials");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Cyclotomic polynomial Φ₁₂(x) = x⁴ - x² + 1
    let cyclotomic_12 = poly_z(&[1, 0, -1, 0, 1]);
    println!("  Cyclotomic Φ₁₂(x) = x⁴ - x² + 1");
    let result = van_hoeij_factor(&cyclotomic_12);
    println!("  Irreducible over Z: {}", result.factors.len() == 1);
    println!();

    // Sophie Germain: x⁴ + 4y⁴ = (x² + 2y² + 2xy)(x² + 2y² - 2xy)
    // Simplified version: x⁴ + 4 (y=1)
    let sophie = poly_z(&[4, 0, 0, 0, 1]);
    println!("  Sophie Germain: x⁴ + 4");
    println!("  Expected: (x² + 2x + 2)(x² - 2x + 2)");
    let result = van_hoeij_factor(&sophie);
    println!("  Found {} factors", result.factors.len());
    println!();

    // Chebyshev polynomial T₄(x) = 8x⁴ - 8x² + 1
    let chebyshev_4 = poly_z(&[1, 0, -8, 0, 8]);
    println!("  Chebyshev T₄(x) = 8x⁴ - 8x² + 1");
    let result = van_hoeij_factor(&chebyshev_4);
    println!("  Irreducible over Z: {}", result.factors.len() == 1);
    println!();

    // Discriminant test polynomial: x³ - 3x + 1 (discriminant = 81)
    let disc_test = poly_z(&[1, -3, 0, 1]);
    println!("  x³ - 3x + 1 (roots related to cos(2π/9))");
    let result = van_hoeij_factor(&disc_test);
    println!("  Irreducible over Z: {}", result.factors.len() == 1);
    println!();
}

/// Example 6: Performance Comparison
fn example_6_performance_comparison() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 6: Performance Benchmarks");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    let iterations = 100;

    // Benchmark (x+1)(x+2)
    let f1 = poly_z(&[2, 3, 1]);
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = van_hoeij_factor(&f1);
    }
    let elapsed = start.elapsed();
    println!("  (x+1)(x+2): {:?} per factorization", elapsed / iterations);

    // Benchmark x⁴ - 1
    let f2 = poly_z(&[-1, 0, 0, 0, 1]);
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = van_hoeij_factor(&f2);
    }
    let elapsed = start.elapsed();
    println!("  x⁴ - 1:    {:?} per factorization", elapsed / iterations);

    // Benchmark bivariate
    let f3 = sparse_poly(&[(1, &[2, 0]), (3, &[1, 0]), (2, &[0, 0])], 2);
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = factor_multivariate(&f3);
    }
    let elapsed = start.elapsed();
    println!("  Bivariate: {:?} per factorization", elapsed / iterations);

    // Benchmark trivariate
    let f4 = sparse_poly(&[(1, &[2, 0, 0]), (3, &[1, 0, 0]), (2, &[0, 0, 0])], 3);
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = factor_multivariate(&f4);
    }
    let elapsed = start.elapsed();
    println!("  Trivariate: {:?} per factorization", elapsed / iterations);
    println!();
}

/// Helper to print factors
fn print_factors(result: &MultivariateFactorResult) {
    for (i, f) in result.factors.iter().enumerate() {
        print!("    Factor {}: ", i + 1);
        let terms: Vec<String> = f
            .terms()
            .iter()
            .map(|(m, c)| {
                let coeff = format!("{:?}", c.0);
                let num_vars = f.num_vars();
                let mut vars = String::new();
                for v in 0..num_vars {
                    let exp = m.exponent(v);
                    if exp > 0 {
                        let var_name = match v {
                            0 => "x",
                            1 => "y",
                            2 => "z",
                            _ => "w",
                        };
                        if exp == 1 {
                            vars.push_str(var_name);
                        } else {
                            vars.push_str(&format!("{}^{}", var_name, exp));
                        }
                    }
                }
                if vars.is_empty() {
                    coeff
                } else if coeff == "1" {
                    vars
                } else if coeff == "-1" {
                    format!("-{}", vars)
                } else {
                    format!("{}*{}", coeff, vars)
                }
            })
            .collect();
        println!("{}", terms.join(" + "));
    }
}
