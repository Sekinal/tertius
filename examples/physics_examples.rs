//! Physics and Math Examples using Tertius
//!
//! Run with: cargo run --example physics_examples

use tertius_factor::van_hoeij_factor;
use tertius_groebner::m5gb::groebner_basis;
use tertius_groebner::monomial::PackedMonomial;
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_rings::finite_field::FiniteField;
use tertius_rings::integers::Z;
use tertius_rings::traits::Field;
use tertius_solve::fglm::fglm_convert;

// Helper to create integers
fn z(n: i64) -> Z {
    Z(Integer::new(n))
}

// Helper to create polynomials over Z
fn poly_z(coeffs: &[i64]) -> DensePoly<Z> {
    DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
}

// Finite field GF(101)
type GF101 = FiniteField<101>;

fn ff(n: i64) -> GF101 {
    GF101::new(n.rem_euclid(101) as u64)
}

fn mono(exps: &[u16]) -> PackedMonomial {
    PackedMonomial::new(exps)
}

fn main() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║          Tertius: Physics & Math Examples                  ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");

    example_1_polynomial_arithmetic();
    example_2_factorization();
    example_3_critical_points();
    example_4_intersection();
    example_5_modular_arithmetic();
}

/// Example 1: Basic Polynomial Arithmetic
fn example_1_polynomial_arithmetic() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 1: Polynomial Arithmetic");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Create polynomials: p(x) = x² + 2x + 1, q(x) = x - 1
    let p = poly_z(&[1, 2, 1]); // 1 + 2x + x²
    let q = poly_z(&[-1, 1]);   // -1 + x

    println!("  p(x) = x² + 2x + 1 = (x + 1)²");
    println!("  q(x) = x - 1\n");

    // Multiply
    let product = p.mul(&q);
    println!("  p(x) × q(x) = {:?}", product.coeffs());
    println!("             = x³ + x² - x - 1\n");

    // Add
    let sum = p.add(&q);
    println!("  p(x) + q(x) = {:?}", sum.coeffs());
    println!("             = x² + 3x\n");
}

/// Example 2: Polynomial Factorization
fn example_2_factorization() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 2: Polynomial Factorization (Van Hoeij + LLL)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Factor x⁴ - 1 = (x-1)(x+1)(x²+1)
    let p1 = poly_z(&[-1, 0, 0, 0, 1]); // x⁴ - 1
    println!("  Factoring: x⁴ - 1");
    let result1 = van_hoeij_factor(&p1);
    println!("  Factors found: {} irreducible factors", result1.factors.len());
    for (i, f) in result1.factors.iter().enumerate() {
        println!("    Factor {}: {:?}", i + 1, f.coeffs());
    }
    println!();

    // Factor x⁴ + 4 (Sophie Germain identity)
    // x⁴ + 4 = (x² + 2x + 2)(x² - 2x + 2)
    let p2 = poly_z(&[4, 0, 0, 0, 1]); // x⁴ + 4
    println!("  Factoring: x⁴ + 4 (Sophie Germain identity)");
    let result2 = van_hoeij_factor(&p2);
    println!("  Factors found: {} irreducible factors", result2.factors.len());
    for (i, f) in result2.factors.iter().enumerate() {
        println!("    Factor {}: {:?}", i + 1, f.coeffs());
    }
    println!();

    // Factor (x+1)⁵
    let p3 = poly_z(&[1, 5, 10, 10, 5, 1]); // (x+1)⁵
    println!("  Factoring: (x+1)⁵ = x⁵ + 5x⁴ + 10x³ + 10x² + 5x + 1");
    let result3 = van_hoeij_factor(&p3);
    println!("  Factors found: {} factors (with multiplicity)", result3.factors.len());
    println!();
}

/// Example 3: Finding Critical Points (Physics: Equilibrium)
fn example_3_critical_points() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 3: Critical Points of a Potential (Equilibrium)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    println!("  Problem: Find critical points of V(x,y) = x² + y² - xy");
    println!("  Setting ∇V = 0:");
    println!("    ∂V/∂x = 2x - y = 0");
    println!("    ∂V/∂y = 2y - x = 0\n");

    // System: 2x - y = 0, 2y - x = 0
    // In GF(101): 2x + 100y = 0, 100x + 2y = 0
    let system = vec![
        // 2x - y (= 2x + 100y in GF(101))
        vec![(ff(2), mono(&[1, 0])), (ff(100), mono(&[0, 1]))],
        // 2y - x (= 100x + 2y in GF(101))
        vec![(ff(100), mono(&[1, 0])), (ff(2), mono(&[0, 1]))],
    ];

    let gb = groebner_basis(system);
    println!("  Gröbner basis computed: {} elements", gb.len());

    let result = fglm_convert(&gb, 2);
    println!("  Lex basis (triangular form): {} polynomials", result.lex_basis.len());
    println!("  Quotient dimension: {} (number of solutions)", result.dimension);
    println!();
    println!("  → Solution: (x, y) = (0, 0) - the unique critical point!\n");
}

/// Example 4: Intersection of Curves
fn example_4_intersection() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 4: Intersection of Circle and Line");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    println!("  Problem: Find intersection of x² + y² = 1 and y = x");
    println!("  Substituting: x² + x² = 1 → 2x² = 1 → x = ±1/√2\n");

    // System: x² + y² - 1 = 0, y - x = 0
    // Over GF(101): x² + y² + 100 = 0, y + 100x = 0
    let system = vec![
        // x² + y² - 1
        vec![
            (ff(1), mono(&[2, 0])),
            (ff(1), mono(&[0, 2])),
            (ff(100), mono(&[0, 0])),
        ],
        // y - x
        vec![(ff(1), mono(&[0, 1])), (ff(100), mono(&[1, 0]))],
    ];

    let gb = groebner_basis(system);
    println!("  Gröbner basis computed: {} elements", gb.len());

    let result = fglm_convert(&gb, 2);
    println!("  Quotient dimension: {} solutions", result.dimension);
    println!();
    println!("  → Two intersection points: (1/√2, 1/√2) and (-1/√2, -1/√2)\n");
}

/// Example 5: Modular Arithmetic (Cryptography/Physics)
fn example_5_modular_arithmetic() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Example 5: Modular Arithmetic in GF(101)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    let a = ff(42);
    let b = ff(17);

    println!("  Working in GF(101) (integers mod 101)\n");

    // Addition
    let sum = a + b;
    println!("  42 + 17 = {} (mod 101)", sum.value());

    // Multiplication
    let prod = a * b;
    println!("  42 × 17 = {} (mod 101)", prod.value());

    // Multiplicative inverse
    let a_inv = a.inv().unwrap();
    println!("  42⁻¹ = {} (mod 101)", a_inv.value());

    // Verify inverse
    let check = a * a_inv;
    println!("  42 × 42⁻¹ = {} (should be 1)", check.value());

    // Division
    let div = a * b.inv().unwrap();
    println!("  42 / 17 = {} (mod 101)", div.value());

    // Powers using repeated squaring
    let power = {
        let mut result = ff(1);
        let mut base = a;
        let mut exp = 10u32;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base;
            }
            base = base * base;
            exp /= 2;
        }
        result
    };
    println!("  42¹⁰ = {} (mod 101)", power.value());

    println!();
    println!("  → Finite field arithmetic: exact, no floating point errors!\n");
}
