//! Integration over exponential extensions.
//!
//! For θ = exp(u), we have θ' = u'·θ.
//!
//! Key property: θ cannot appear in denominators of integrals
//! (the "exponential Liouville theorem").
//!
//! # Algorithm
//!
//! To integrate f(θ) in K(θ) where θ = exp(u):
//!
//! 1. Write f = p/q where p, q ∈ K[θ]
//! 2. If deg(q) > 0 in θ, need special handling (residues)
//! 3. For polynomials: match coefficients with θ' = u'θ

use tertius_poly::dense::DensePoly;
use tertius_rings::traits::Field;

/// Result of integrating over an exponential extension.
#[derive(Clone, Debug)]
pub enum ExpExtIntegrationResult<F: Field> {
    /// Found an antiderivative in K(θ)
    Success {
        /// Coefficients of θⁱ in the result (may include negative powers)
        coeffs: Vec<(i32, F)>, // (power of θ, coefficient from K)
    },
    /// The integrand has no elementary antiderivative
    NonElementary,
    /// Algorithm couldn't complete
    Failed,
}

/// Integrates a polynomial in θ where θ = exp(u) and θ' = u'·θ.
///
/// For f = aₙθⁿ + ... + a₁θ + a₀ ∈ K[θ], we seek g ∈ K[θ] with g' = f.
///
/// Key observation: d/dx(bθⁿ) = b'θⁿ + b·n·u'·θⁿ = (b' + n·u'·b)θⁿ
///
/// So to get aₙθⁿ, we need: b' + n·u'·b = aₙ
///
/// This is a first-order linear ODE in b!
pub fn integrate_exp_polynomial<F: Field>(
    f: &DensePoly<F>,
    u_deriv: &F, // This is u' (the derivative of the exponent)
) -> ExpExtIntegrationResult<F> {
    if f.is_zero() {
        return ExpExtIntegrationResult::Success { coeffs: vec![] };
    }

    let n = f.degree();
    let mut result_coeffs = Vec::new();

    // For each term aᵢθⁱ in f, we need to solve: b' + i·u'·b = aᵢ
    //
    // If aᵢ and u' are both constants (in Q), then:
    // - For i = 0: b' = a₀ → b = a₀·x (requires derivative of constants = 0)
    // - For i ≠ 0: constant b works if i·u'·b = aᵢ → b = aᵢ/(i·u')
    //   But b' = 0, so we need a₀ = 0 for this to work OR proper integration

    // Simplified case: u' is a non-zero constant
    // Then for θⁱ term with i ≠ 0:
    // d/dx(b·θⁱ) = b'·θⁱ + b·i·u'·θⁱ = (b' + i·u'·b)·θⁱ
    // If b is constant: b' = 0, so we need b = aᵢ/(i·u')

    let u_prime = u_deriv.clone();

    if u_prime == F::zero() {
        // u' = 0 means θ = exp(constant), which is itself constant
        // Then we're just integrating a polynomial in a constant θ
        // This degenerates to the base case
        return ExpExtIntegrationResult::Failed;
    }

    for i in 0..=n {
        let ai = f.coeff(i);

        if ai == F::zero() {
            continue;
        }

        if i == 0 {
            // Constant term: need ∫a₀ dx
            // This requires integration in the base field
            // For a₀ constant, we'd get a₀·x, but x is not in K(θ)
            // This means we need to handle it separately
            // For now, return Failed for non-zero constant terms
            return ExpExtIntegrationResult::Failed;
        }

        // For i ≠ 0: bᵢ = aᵢ/(i·u')
        let i_elem = F::one().mul_by_scalar(i as i64);
        let factor = i_elem * u_prime.clone();
        let bi = ai * factor.inv().unwrap();

        result_coeffs.push((i as i32, bi));
    }

    ExpExtIntegrationResult::Success { coeffs: result_coeffs }
}

/// Integrates θⁿ where θ = exp(x).
///
/// ∫ exp(nx) dx = exp(nx)/n for n ≠ 0
/// ∫ 1 dx = x
pub fn integrate_exp_power(n: i32) -> ExpPowerIntegral {
    if n == 0 {
        ExpPowerIntegral::Linear // Result is x
    } else {
        ExpPowerIntegral::ExpTerm {
            exp_coeff: n,
            scalar: (1, n), // Result is exp(nx)/n
        }
    }
}

/// The integral of θⁿ where θ = exp(x).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExpPowerIntegral {
    /// The result is x (for n = 0)
    Linear,
    /// The result is (num/denom) * exp(exp_coeff * x)
    ExpTerm {
        exp_coeff: i32,
        scalar: (i64, i32), // numerator, denominator
    },
}

/// Special integration formulas for products of polynomials and exponentials.
///
/// ∫ xⁿ exp(ax) dx can be computed using integration by parts.
/// Result: exp(ax) * Σᵢ₌₀ⁿ (-1)ⁱ n!/((n-i)! aⁱ⁺¹) xⁿ⁻ⁱ
pub fn integrate_poly_times_exp(poly_deg: usize, exp_coeff: i32) -> Vec<(i64, usize, i64)> {
    // Returns terms of form (coeff, power of x, power of a in denominator)
    // Result: exp(ax) * Σ terms

    if exp_coeff == 0 {
        // Degenerate case: just integrate xⁿ
        // ∫ xⁿ dx = xⁿ⁺¹/(n+1)
        return vec![(1, poly_deg + 1, 0)];
    }

    let n = poly_deg;
    let a = exp_coeff;
    let mut terms = Vec::new();

    // Using the formula from integration by parts:
    // ∫ xⁿ eᵃˣ dx = eᵃˣ Σᵢ₌₀ⁿ (-1)ⁱ n!/((n-i)! aⁱ⁺¹) xⁿ⁻ⁱ

    let mut n_factorial: i64 = 1;
    for k in 1..=n {
        n_factorial *= k as i64;
    }

    for i in 0..=n {
        // Term: (-1)^i * n! / ((n-i)! * a^(i+1)) * x^(n-i)
        let sign: i64 = if i % 2 == 0 { 1 } else { -1 };

        let mut n_minus_i_factorial: i64 = 1;
        for k in 1..=(n - i) {
            n_minus_i_factorial *= k as i64;
        }

        let coeff = sign * n_factorial / n_minus_i_factorial;
        let x_power = n - i;
        let a_power = (i + 1) as i64;

        // We store power of a in denominator as positive
        // The actual denominator is a^(i+1)
        terms.push((coeff, x_power, a_power));
    }

    terms
}

/// Checks if a rational function in θ = exp(u) is integrable.
///
/// Key theorem: If θ = exp(u), then ∫(p/q) elementary implies q | p
/// in a specific sense (the denominator cannot grow).
pub fn is_exp_extension_integrable<F: Field>(
    numerator: &DensePoly<F>,
    denominator: &DensePoly<F>,
) -> Option<bool> {
    // If denominator has positive degree in θ, check special conditions
    if denominator.degree() > 0 {
        // In general, rational functions in exp with non-trivial denominator
        // require careful analysis
        return None;
    }

    // Polynomial in θ: analyze each term
    if numerator.is_zero() {
        return Some(true);
    }

    // For polynomials in exp(u), integrability depends on u'
    // We'd need to check if the first-order ODEs are solvable
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn poly(coeffs: &[i64]) -> DensePoly<Q> {
        DensePoly::new(coeffs.iter().map(|&n| q(n)).collect())
    }

    #[test]
    fn test_integrate_zero() {
        let zero = DensePoly::zero();
        let u_deriv = q(1);

        let result = integrate_exp_polynomial(&zero, &u_deriv);

        if let ExpExtIntegrationResult::Success { coeffs } = result {
            assert!(coeffs.is_empty());
        } else {
            panic!("Expected success");
        }
    }

    #[test]
    fn test_integrate_exp_x() {
        // f = θ = exp(x), so coeffs = [0, 1]
        // u = x, u' = 1
        // ∫ exp(x) dx = exp(x), so result should be θ with coefficient 1

        let f = poly(&[0, 1]); // θ¹
        let u_deriv = q(1);

        let result = integrate_exp_polynomial(&f, &u_deriv);

        if let ExpExtIntegrationResult::Success { coeffs } = result {
            assert_eq!(coeffs.len(), 1);
            assert_eq!(coeffs[0], (1, q(1))); // 1 * θ¹
        } else {
            panic!("Expected success");
        }
    }

    #[test]
    fn test_integrate_exp_2x() {
        // f = θ² where θ = exp(x)
        // This represents exp(2x) = (exp(x))²
        // ∫ exp(2x) dx = exp(2x)/2, so result is θ²/2

        let f = poly(&[0, 0, 1]); // θ²
        let u_deriv = q(1); // u = x, u' = 1

        let result = integrate_exp_polynomial(&f, &u_deriv);

        if let ExpExtIntegrationResult::Success { coeffs } = result {
            assert_eq!(coeffs.len(), 1);
            assert_eq!(coeffs[0], (2, q(1) * q(2).inv().unwrap())); // (1/2) * θ²
        } else {
            panic!("Expected success");
        }
    }

    #[test]
    fn test_exp_power_integral() {
        assert_eq!(integrate_exp_power(0), ExpPowerIntegral::Linear);

        let result = integrate_exp_power(2);
        assert_eq!(result, ExpPowerIntegral::ExpTerm {
            exp_coeff: 2,
            scalar: (1, 2),
        });
    }

    #[test]
    fn test_poly_times_exp() {
        // ∫ x exp(x) dx = exp(x)(x - 1)
        let terms = integrate_poly_times_exp(1, 1);

        // Should have two terms: coeff*x^1/a^1 and coeff*x^0/a^2
        assert_eq!(terms.len(), 2);
        // x term: (+1) * 1!/0! * x¹ / a¹ = x/a
        assert_eq!(terms[0], (1, 1, 1));
        // constant term: (-1) * 1!/1! * x⁰ / a² = -1/a²
        assert_eq!(terms[1], (-1, 0, 2));
    }

    #[test]
    fn test_poly_times_exp_squared() {
        // ∫ x² exp(x) dx = exp(x)(x² - 2x + 2)
        let terms = integrate_poly_times_exp(2, 1);

        assert_eq!(terms.len(), 3);
        // x² term: 2!/2! * x² / a¹ = x²
        assert_eq!(terms[0], (1, 2, 1));
        // x term: (-1) * 2!/1! / a² = -2x/a²
        assert_eq!(terms[1], (-2, 1, 2));
        // constant: 2!/0! / a³ = 2/a³
        assert_eq!(terms[2], (2, 0, 3));
    }
}
