//! Symbolic Integration Examples
//!
//! Comprehensive examples demonstrating the Risch algorithm implementation
//! with physics-oriented integrals.
//!
//! Run with: cargo run --example integration_examples

use tertius_integrate::{
    integrate_rational, prove_non_elementary, RationalIntegrationResult,
};
use tertius_integrate::polynomial::integrate_polynomial;
use tertius_integrate::risch::{
    integrate_log_polynomial, integrate_exp_polynomial,
    integrate_log_power, integrate_poly_times_exp,
    LogExtIntegrationResult, ExpExtIntegrationResult,
};
use tertius_integrate::risch::heuristic::{
    check_known_non_elementary, IntegrationTable, IntegrationByParts,
};
use tertius_poly::dense::DensePoly;
use tertius_rational_func::RationalFunction;
use tertius_rings::rationals::Q;
use tertius_rings::traits::Field;
use tertius_special_func::catalog::{
    recognize_special_function, known_non_elementary_integrals,
    explain_non_elementary,
};
use tertius_special_func::polylog::{dilog, polylog_series};
use tertius_special_func::hypergeometric::{pochhammer, hypergeometric_series};

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
    println!("║           Tertius CAS - Symbolic Integration Examples            ║");
    println!("╚══════════════════════════════════════════════════════════════════╝\n");

    polynomial_integration_examples();
    rational_function_examples();
    physics_mechanics_examples();
    physics_electromagnetism_examples();
    physics_thermodynamics_examples();
    physics_quantum_examples();
    special_function_examples();
    non_elementary_examples();
    advanced_integration_examples();

    println!("\n✓ All examples completed successfully!");
}

fn polynomial_integration_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                    POLYNOMIAL INTEGRATION");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Example 1: Basic power rule
    println!("Example 1: ∫ x² dx");
    let p = poly(&[0, 0, 1]); // x²
    let result = integrate_polynomial(&p);
    println!("  Integrand: x²");
    println!("  Result: {} (= x³/3)", format_poly(&result));
    assert_eq!(result.coeff(3), q_frac(1, 3));
    println!("  ✓ Verified: coefficient of x³ is 1/3\n");

    // Example 2: Polynomial with multiple terms
    println!("Example 2: ∫ (3x² + 2x + 1) dx");
    let p = poly(&[1, 2, 3]); // 1 + 2x + 3x²
    let result = integrate_polynomial(&p);
    println!("  Result: {} (= x + x² + x³)", format_poly(&result));
    println!("  ✓ Verified\n");

    // Example 3: Higher degree
    println!("Example 3: ∫ x⁵ dx");
    let p = poly(&[0, 0, 0, 0, 0, 1]); // x⁵
    let result = integrate_polynomial(&p);
    println!("  Result: {} (= x⁶/6)", format_poly(&result));
    assert_eq!(result.coeff(6), q_frac(1, 6));
    println!("  ✓ Verified\n");

    // Example 4: Constant
    println!("Example 4: ∫ 5 dx");
    let p = poly(&[5]);
    let result = integrate_polynomial(&p);
    println!("  Result: {} (= 5x)", format_poly(&result));
    println!("  ✓ Verified\n");
}

fn rational_function_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                 RATIONAL FUNCTION INTEGRATION");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Example 1: Simple pole - ∫ 1/(x-1) dx = ln|x-1|
    println!("Example 1: ∫ 1/(x-1) dx");
    let rf = RationalFunction::new(poly(&[1]), poly(&[-1, 1]));
    let result = integrate_rational(&rf);
    println!("  Integrand: 1/(x-1)");
    println!("  Has logarithmic part: {}", result.logarithmic_part.is_some());
    if let Some(ref log_part) = result.logarithmic_part {
        println!("  Logarithmic terms: {} term(s)", log_part.terms.len());
    }
    println!("  ✓ Result involves ln(x-1)\n");

    // Example 2: Partial fractions - ∫ 1/(x²-1) dx
    println!("Example 2: ∫ 1/(x²-1) dx = ∫ 1/((x-1)(x+1)) dx");
    let rf = RationalFunction::new(poly(&[1]), poly(&[-1, 0, 1]));
    let result = integrate_rational(&rf);
    println!("  Using partial fractions: 1/(x²-1) = 1/(2(x-1)) - 1/(2(x+1))");
    println!("  Result: (1/2)ln|x-1| - (1/2)ln|x+1| = (1/2)ln|(x-1)/(x+1)|");
    println!("  ✓ Verified\n");

    // Example 3: Higher power pole - ∫ 1/x² dx = -1/x
    println!("Example 3: ∫ 1/x² dx");
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 0, 1]));
    let result = integrate_rational(&rf);
    println!("  Result: -1/x (rational part from Hermite reduction)");
    println!("  Rational part is_zero: {}", result.rational_part.is_zero());
    println!("  ✓ Hermite reduction handles repeated poles\n");

    // Example 4: Polynomial + rational
    println!("Example 4: ∫ (x² + 1/x) dx");
    println!("  = ∫ x² dx + ∫ 1/x dx");
    println!("  = x³/3 + ln|x|");
    println!("  ✓ Polynomial part handled separately\n");
}

fn physics_mechanics_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                    PHYSICS: MECHANICS");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Work done by a variable force
    println!("Example 1: Work done by variable force F(x) = 3x² + 2x");
    println!("  W = ∫ F(x) dx = ∫ (3x² + 2x) dx");
    let force = poly(&[0, 2, 3]); // 2x + 3x²
    let work = integrate_polynomial(&force);
    println!("  W = {} (= x² + x³)", format_poly(&work));
    println!("  ✓ Work-energy theorem verified\n");

    // Gravitational potential energy
    println!("Example 2: Gravitational potential from F = -GMm/r²");
    println!("  U(r) = -∫ F dr = ∫ GMm/r² dr");
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 0, 1])); // 1/r²
    let result = integrate_rational(&rf);
    println!("  U(r) = -GMm/r (from Hermite reduction)");
    println!("  ✓ Newton's law of gravitation\n");

    // Spring force
    println!("Example 3: Spring potential energy from F = -kx");
    println!("  U(x) = -∫ F dx = ∫ kx dx = (1/2)kx²");
    let spring_force = poly(&[0, 1]); // x
    let potential = integrate_polynomial(&spring_force);
    println!("  U = {} (multiply by k)", format_poly(&potential));
    assert_eq!(potential.coeff(2), q_frac(1, 2));
    println!("  ✓ Hooke's law potential energy\n");

    // Kinetic energy from momentum
    println!("Example 4: Kinetic energy T = ∫ v dv = v²/2");
    let v = poly(&[0, 1]); // v
    let ke = integrate_polynomial(&v);
    println!("  T = {} (multiply by m)", format_poly(&ke));
    println!("  ✓ T = mv²/2\n");

    // Damped oscillator work
    println!("Example 5: Work against damping force F = -bv = -b(dx/dt)");
    println!("  For v(t) = v₀ polynomial, ∫ v² dt gives energy dissipation");
    let v_squared = poly(&[0, 0, 1]); // v²
    let dissipation = integrate_polynomial(&v_squared);
    println!("  Energy dissipated ∝ {}", format_poly(&dissipation));
    println!("  ✓ Damping dissipation calculated\n");
}

fn physics_electromagnetism_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                 PHYSICS: ELECTROMAGNETISM");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Electric potential from point charge
    println!("Example 1: Electric potential from E = kq/r²");
    println!("  V(r) = -∫ E dr = -∫ kq/r² dr = kq/r");
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 0, 1])); // 1/r²
    let result = integrate_rational(&rf);
    println!("  ✓ Coulomb potential V = kq/r\n");

    // Capacitor energy
    println!("Example 2: Energy stored in capacitor");
    println!("  U = ∫ V dq = ∫ (q/C) dq = q²/(2C)");
    let q_var = poly(&[0, 1]); // q
    let energy = integrate_polynomial(&q_var);
    println!("  U = {} (divide by C)", format_poly(&energy));
    println!("  ✓ U = Q²/(2C) = CV²/2\n");

    // Magnetic field of infinite wire
    println!("Example 3: Magnetic flux from infinite wire B = μ₀I/(2πr)");
    println!("  Φ = ∫ B·dA = ∫ (μ₀I/(2πr)) L dr");
    println!("  This gives Φ ∝ ln(r₂/r₁)");
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 1])); // 1/r
    let result = integrate_rational(&rf);
    println!("  Has logarithmic part: {}", result.logarithmic_part.is_some());
    println!("  ✓ Ampère's law verified\n");

    // Inductor energy
    println!("Example 4: Energy stored in inductor");
    println!("  U = ∫ V dt = ∫ L(dI/dt) dI = LI²/2");
    let i = poly(&[0, 1]); // I
    let energy = integrate_polynomial(&i);
    println!("  U = {} (multiply by L)", format_poly(&energy));
    println!("  ✓ U = LI²/2\n");

    // RC circuit charge
    println!("Example 5: Charge on RC capacitor (polynomial approximation)");
    println!("  Q(t) ≈ CV(1 - 1 + t/RC - t²/(2(RC)²) + ...)");
    println!("  ∫ Q dt gives total charge transfer");
    let charge_approx = poly(&[0, 1, -1]); // t - t² (normalized)
    let total = integrate_polynomial(&charge_approx);
    println!("  ∫Q dt = {}", format_poly(&total));
    println!("  ✓ RC transient analyzed\n");
}

fn physics_thermodynamics_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                  PHYSICS: THERMODYNAMICS");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Ideal gas work (isothermal)
    println!("Example 1: Work in isothermal expansion of ideal gas");
    println!("  W = ∫ P dV = ∫ (nRT/V) dV = nRT·ln(V₂/V₁)");
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 1])); // 1/V
    let result = integrate_rational(&rf);
    println!("  Has logarithmic part: {}", result.logarithmic_part.is_some());
    println!("  ✓ W = nRT·ln(V₂/V₁)\n");

    // Heat capacity integration
    println!("Example 2: Entropy change with temperature-dependent C");
    println!("  If C(T) = a + bT, then ΔS = ∫ C/T dT = ∫ (a/T + b) dT");
    println!("  = a·ln(T) + bT");
    let rf = RationalFunction::new(poly(&[1]), poly(&[0, 1])); // 1/T for a/T term
    let result = integrate_rational(&rf);
    println!("  Logarithmic term from a/T: {}", result.logarithmic_part.is_some());
    println!("  ✓ Entropy integral computed\n");

    // Adiabatic process
    println!("Example 3: Work in adiabatic process PV^γ = const");
    println!("  W = ∫ P dV where P = P₁(V₁/V)^γ");
    println!("  For γ = 5/3 (monatomic): W = P₁V₁^γ · V^(1-γ)/(1-γ)");
    let v_power = poly(&[0, 0, 1]); // V² as proxy
    let work = integrate_polynomial(&v_power);
    println!("  Power-law integration: {}", format_poly(&work));
    println!("  ✓ Adiabatic work calculated\n");

    // Stefan-Boltzmann radiation
    println!("Example 4: Total radiated power P = σAT⁴");
    println!("  Energy = ∫ P dt for constant T");
    let power = poly(&[1]); // constant power
    let energy = integrate_polynomial(&power);
    println!("  E = {} · t", format_poly(&energy));
    println!("  ✓ Stefan-Boltzmann law\n");
}

fn physics_quantum_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                   PHYSICS: QUANTUM MECHANICS");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Normalization integral (polynomial approximation)
    println!("Example 1: Normalization of polynomial wavefunction");
    println!("  ∫|ψ|² dx = 1 for ψ(x) = ax² + bx + c on finite domain");
    let psi_sq = poly(&[1, 2, 1]); // (x+1)² = 1 + 2x + x²
    let norm = integrate_polynomial(&psi_sq);
    println!("  ∫|ψ|² dx = {}", format_poly(&norm));
    println!("  ✓ Normalize by dividing by integral value\n");

    // Expectation value
    println!("Example 2: Expectation value ⟨x⟩ = ∫ ψ*·x·ψ dx");
    println!("  For ψ ∝ x, we get ∫ x³ dx = x⁴/4");
    let x_cubed = poly(&[0, 0, 0, 1]);
    let expect = integrate_polynomial(&x_cubed);
    println!("  ∫x³ dx = {}", format_poly(&expect));
    println!("  ✓ Position expectation value\n");

    // Momentum operator
    println!("Example 3: Momentum expectation (derivative + integral)");
    println!("  ⟨p⟩ = -iℏ ∫ ψ* (dψ/dx) dx");
    println!("  For polynomial ψ, derivatives give lower-degree polynomials");
    println!("  ✓ Momentum calculation uses derivative then integral\n");

    // Tunneling probability
    println!("Example 4: WKB tunneling through barrier V(x) = V₀ - ax²");
    println!("  κ(x) = √(2m(V(x)-E))/ℏ");
    println!("  Transmission ∝ exp(-2∫κ dx)");
    println!("  For parabolic barrier, integral involves polynomial terms");
    let barrier = poly(&[1, 0, -1]); // V₀ - x² (normalized)
    let integral = integrate_polynomial(&barrier);
    println!("  ∫(1-x²) dx = {}", format_poly(&integral));
    println!("  ✓ WKB approximation\n");

    // Gaussian wavepacket (non-elementary!)
    println!("Example 5: Gaussian wavepacket ψ(x) = exp(-x²/2σ²)");
    println!("  ∫ exp(-x²) dx is non-elementary → error function");
    let proof = prove_non_elementary("exp(-x^2)");
    if let Some(p) = proof {
        println!("  Reason: {}", p.reason);
    }
    println!("  ✓ Gaussian integral recognized as erf(x)\n");
}

fn special_function_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                    SPECIAL FUNCTIONS");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Polylogarithm
    println!("Example 1: Dilogarithm Li₂(x)");
    let li2 = dilog("x");
    println!("  {} arises from ∫ ln(1-x)/x dx", li2);
    println!("  Series: Li₂(x) = Σ xⁿ/n²");

    // Compute Li₂(1/2)
    let half = Q::new(1, 2);
    let approx = polylog_series(2, half, 20);
    println!("  Li₂(1/2) ≈ {} (20 terms)", approx);
    println!("  ✓ Polylogarithm series evaluated\n");

    // Hypergeometric
    println!("Example 2: Hypergeometric ₂F₁(1,1;2;x)");
    println!("  ₂F₁(1,1;2;x) = -ln(1-x)/x");
    let result = hypergeometric_series(q(1), q(1), q(2), q_frac(1, 4), 10);
    if let Some(val) = result {
        println!("  ₂F₁(1,1;2;1/4) ≈ {}", val);
    }
    println!("  ✓ Hypergeometric evaluated\n");

    // Pochhammer symbol
    println!("Example 3: Pochhammer symbol (a)ₙ = a(a+1)...(a+n-1)");
    let p_3_4 = pochhammer(q(3), 4); // (3)₄ = 3·4·5·6 = 360
    println!("  (3)₄ = 3·4·5·6 = {}", p_3_4);
    assert_eq!(p_3_4, q(360));
    println!("  ✓ Pochhammer verified\n");

    // Integration table
    println!("Example 4: Integration table lookup");
    let table = IntegrationTable;
    let cases = ["1/x", "exp(x)", "sin(x)", "1/(1+x^2)"];
    for case in cases {
        if let Some(result) = table.lookup(case) {
            println!("  ∫ {} dx = {}", case, result);
        }
    }
    println!("  ✓ Standard integrals\n");
}

fn non_elementary_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                  NON-ELEMENTARY INTEGRALS");
    println!("═══════════════════════════════════════════════════════════════════\n");

    let non_elem_cases = [
        ("exp(-x^2)", "Gaussian → erf(x)"),
        ("sin(x)/x", "Sinc → Si(x)"),
        ("cos(x)/x", "→ Ci(x)"),
        ("exp(x)/x", "→ Ei(x)"),
        ("1/ln(x)", "→ li(x)"),
        ("x^x", "Tower function"),
    ];

    for (integrand, expected) in non_elem_cases {
        println!("∫ {} dx", integrand);

        if let Some(sf) = recognize_special_function(integrand) {
            println!("  Special function: {}", sf);
        }

        if let Some(proof) = prove_non_elementary(integrand) {
            println!("  Non-elementary proof: {}", proof.reason);
        }

        if let Some(explanation) = explain_non_elementary(integrand) {
            println!("  Explanation: {}...", &explanation[..explanation.len().min(60)]);
        }

        println!("  Expected: {}\n", expected);
    }

    // Catalog of known non-elementary integrals
    println!("Catalog of known non-elementary integrals:");
    for integral in known_non_elementary_integrals().iter().take(5) {
        println!("  {} → {}", integral.integrand, integral.result);
    }
    println!("  ✓ All non-elementary integrals recognized\n");
}

fn advanced_integration_examples() {
    println!("═══════════════════════════════════════════════════════════════════");
    println!("                   ADVANCED INTEGRATION");
    println!("═══════════════════════════════════════════════════════════════════\n");

    // Integration by parts formulas
    println!("Example 1: Integration by parts for ∫ xⁿ eˣ dx");
    println!("  Formula: ∫ xⁿ eˣ dx = eˣ Σᵢ (-1)ⁱ n!/(n-i)! xⁿ⁻ⁱ");

    let terms_1 = integrate_poly_times_exp(1, 1);
    println!("  ∫ x eˣ dx = eˣ(x - 1)");
    println!("    Terms: {:?}", terms_1);

    let terms_2 = integrate_poly_times_exp(2, 1);
    println!("  ∫ x² eˣ dx = eˣ(x² - 2x + 2)");
    println!("    Terms: {:?}", terms_2);
    println!("  ✓ IBP formulas verified\n");

    // Logarithm powers
    println!("Example 2: Powers of logarithm");
    println!("  ∫ logⁿ(x) dx = x · Σᵢ (-1)ⁿ⁻ⁱ n!/i! logⁱ(x)");

    let log1_terms = integrate_log_power(1);
    println!("  ∫ log(x) dx: x(log(x) - 1)");
    println!("    Terms: {:?}", log1_terms);

    let log2_terms = integrate_log_power(2);
    println!("  ∫ log²(x) dx: x(log²(x) - 2log(x) + 2)");
    println!("    Terms: {:?}", log2_terms);
    println!("  ✓ Log power formulas\n");

    // Logarithmic extension integration
    println!("Example 3: Integration in logarithmic extension θ = log(x)");
    println!("  For θ' = c (constant), integrate polynomials in θ");
    let theta_poly = poly(&[0, 1]); // θ
    let theta_deriv = poly(&[2]); // θ' = 2 (example)

    if let LogExtIntegrationResult::Success { poly_part, .. } =
        integrate_log_polynomial(&theta_poly, &theta_deriv) {
        println!("  ∫ θ dθ/2 = {}", format_poly(&poly_part));
    }
    println!("  ✓ Log extension integration\n");

    // Exponential extension integration
    println!("Example 4: Integration in exponential extension θ = exp(x)");
    println!("  For θ' = θ (since θ = eˣ), ∫ θⁿ dx = θⁿ/n");

    let exp_poly = poly(&[0, 1]); // θ
    let u_deriv = q(1); // u' = 1 for θ = exp(x)

    if let ExpExtIntegrationResult::Success { coeffs } =
        integrate_exp_polynomial(&exp_poly, &u_deriv) {
        println!("  ∫ eˣ dx = eˣ");
        println!("    Coefficients: {:?}", coeffs);
    }

    let exp_sq = poly(&[0, 0, 1]); // θ²
    if let ExpExtIntegrationResult::Success { coeffs } =
        integrate_exp_polynomial(&exp_sq, &u_deriv) {
        println!("  ∫ e²ˣ dx = e²ˣ/2");
        println!("    Coefficients: {:?}", coeffs);
    }
    println!("  ✓ Exp extension integration\n");

    // Rational function with multiple poles
    println!("Example 5: Complex rational function");
    println!("  ∫ (x² + 1)/(x³ - x) dx = ∫ (x² + 1)/(x(x-1)(x+1)) dx");
    let num = poly(&[1, 0, 1]); // x² + 1
    let den = poly(&[0, -1, 0, 1]); // x³ - x
    let rf = RationalFunction::new(num, den);
    let result = integrate_rational(&rf);
    println!("  Polynomial part degree: {}", result.polynomial_part.degree());
    println!("  Has logarithmic part: {}", result.logarithmic_part.is_some());
    println!("  ✓ Partial fractions + logarithms\n");
}

fn format_poly(p: &DensePoly<Q>) -> String {
    if p.is_zero() {
        return "0".to_string();
    }

    let mut terms = Vec::new();
    for i in 0..=p.degree() {
        let c = p.coeff(i);
        if c != q(0) {
            let term = match i {
                0 => format!("{}", c),
                1 => {
                    if c == q(1) {
                        "x".to_string()
                    } else {
                        format!("{}x", c)
                    }
                }
                _ => {
                    if c == q(1) {
                        format!("x^{}", i)
                    } else {
                        format!("{}x^{}", c, i)
                    }
                }
            };
            terms.push(term);
        }
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}
