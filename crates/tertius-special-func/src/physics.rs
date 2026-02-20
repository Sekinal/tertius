//! Physics-oriented special functions (starter set).

use std::f64::consts::PI;

/// Legendre polynomial P_n(x).
#[must_use]
pub fn legendre_p(n: u32, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut p_nm2 = 1.0;
            let mut p_nm1 = x;
            for k in 2..=n {
                let kf = f64::from(k);
                let p_n = ((2.0 * kf - 1.0) * x * p_nm1 - (kf - 1.0) * p_nm2) / kf;
                p_nm2 = p_nm1;
                p_nm1 = p_n;
            }
            p_nm1
        }
    }
}

/// Hermite polynomial H_n(x) (physicists' convention).
#[must_use]
pub fn hermite_h(n: u32, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut h_nm2 = 1.0;
            let mut h_nm1 = 2.0 * x;
            for k in 2..=n {
                let h_n = 2.0 * x * h_nm1 - 2.0 * (f64::from(k) - 1.0) * h_nm2;
                h_nm2 = h_nm1;
                h_nm1 = h_n;
            }
            h_nm1
        }
    }
}

/// Laguerre polynomial L_n(x).
#[must_use]
pub fn laguerre_l(n: u32, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 1.0 - x,
        _ => {
            let mut l_nm2 = 1.0;
            let mut l_nm1 = 1.0 - x;
            for k in 2..=n {
                let kf = f64::from(k);
                let l_n = ((2.0 * kf - 1.0 - x) * l_nm1 - (kf - 1.0) * l_nm2) / kf;
                l_nm2 = l_nm1;
                l_nm1 = l_n;
            }
            l_nm1
        }
    }
}

/// Bessel function J_n(x) via truncated power series.
#[must_use]
pub fn bessel_j(n: u32, x: f64, terms: usize) -> f64 {
    let mut sum = 0.0;
    let x_half = x / 2.0;
    for m in 0..terms {
        let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
        let m_u = m as u32;
        let denom = factorial(m_u) * factorial(m_u + n);
        let power = 2 * m as i32 + n as i32;
        let term = sign * x_half.powi(power) / denom;
        sum += term;
    }
    sum
}

/// Real-valued spherical harmonic Y_0^0.
#[must_use]
pub fn spherical_harmonic_y_00() -> f64 {
    1.0 / (2.0 * PI.sqrt())
}

fn factorial(n: u32) -> f64 {
    (1..=n).fold(1.0, |acc, k| acc * f64::from(k))
}

/// Symbolic physics-special-function expression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PhysicsSpecialExpr {
    /// Legendre polynomial `P_n(arg)`.
    Legendre { n: u32, arg: String },
    /// Hermite polynomial `H_n(arg)` (physicists' convention).
    Hermite { n: u32, arg: String },
    /// Laguerre polynomial `L_n(arg)`.
    Laguerre { n: u32, arg: String },
    /// Bessel function `J_n(arg)`.
    BesselJ { n: i32, arg: String },
}

/// Derivative transformation rule as a symbolic string.
#[must_use]
pub fn derivative_rule(_expr: &PhysicsSpecialExpr) -> String {
    match _expr {
        PhysicsSpecialExpr::Legendre { n, arg } => {
            if *n == 0 {
                "0".to_string()
            } else {
                format!(
                    "({n}/({arg}^2 - 1)) * ({arg}*P_{n}({arg}) - P_{}({arg}))",
                    n - 1
                )
            }
        }
        PhysicsSpecialExpr::Hermite { n, arg } => {
            if *n == 0 {
                "0".to_string()
            } else {
                format!("{}*H_{}({arg})", 2 * n, n - 1)
            }
        }
        PhysicsSpecialExpr::Laguerre { n, arg } => {
            if *n == 0 {
                "0".to_string()
            } else {
                format!("-L_{}^{{(1)}}({arg})", n - 1)
            }
        }
        PhysicsSpecialExpr::BesselJ { n, arg } => {
            format!("(J_{}({arg}) - J_{}({arg}))/2", n - 1, n + 1)
        }
    }
}

/// Integral transformation rule when available.
#[must_use]
pub fn integral_rule(_expr: &PhysicsSpecialExpr) -> Option<String> {
    match _expr {
        PhysicsSpecialExpr::Legendre { n, arg } => {
            if *n == 0 {
                Some(format!("{arg} + C"))
            } else {
                Some(format!(
                    "(P_{}({arg}) - P_{}({arg}))/{} + C",
                    n + 1,
                    n - 1,
                    2 * n + 1
                ))
            }
        }
        PhysicsSpecialExpr::BesselJ { n, arg } if *n == 1 => Some(format!("-J_0({arg}) + C")),
        _ => None,
    }
}

/// Simplification rule for special values (`x=0`, `x=±1`, etc.).
#[must_use]
pub fn simplify_rule(_expr: &PhysicsSpecialExpr) -> Option<String> {
    match _expr {
        PhysicsSpecialExpr::Legendre { n: _, arg } if arg == "1" => Some("1".to_string()),
        PhysicsSpecialExpr::Legendre { n, arg } if arg == "-1" => {
            if n % 2 == 0 {
                Some("1".to_string())
            } else {
                Some("-1".to_string())
            }
        }
        PhysicsSpecialExpr::Laguerre { arg, .. } if arg == "0" => Some("1".to_string()),
        PhysicsSpecialExpr::Hermite { n, arg } if arg == "0" => {
            if n % 2 == 1 {
                Some("0".to_string())
            } else {
                let m = n / 2;
                let value = (-1.0f64).powi(m as i32) * factorial(*n) / factorial(m);
                Some(format!("{}", value.round() as i64))
            }
        }
        PhysicsSpecialExpr::BesselJ { n, arg } if arg == "0" => {
            if *n == 0 {
                Some("1".to_string())
            } else {
                Some("0".to_string())
            }
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_legendre_p_low_orders() {
        assert!((legendre_p(0, 0.3) - 1.0).abs() < 1e-12);
        assert!((legendre_p(1, 0.3) - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_h_low_orders() {
        assert!((hermite_h(0, 1.2) - 1.0).abs() < 1e-12);
        assert!((hermite_h(1, 1.2) - 2.4).abs() < 1e-12);
    }

    #[test]
    fn test_laguerre_l_low_orders() {
        assert!((laguerre_l(0, 0.7) - 1.0).abs() < 1e-12);
        assert!((laguerre_l(1, 0.7) - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_bessel_j0_small_x() {
        // J0(0.2) ~ 0.990024972
        let j0 = bessel_j(0, 0.2, 12);
        assert!((j0 - 0.990024972).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_harmonic_y00() {
        let y00 = spherical_harmonic_y_00();
        let expected = 1.0 / (2.0 * std::f64::consts::PI.sqrt());
        assert!((y00 - expected).abs() < 1e-12);
    }

    #[test]
    fn test_derivative_rule_legendre() {
        let expr = PhysicsSpecialExpr::Legendre {
            n: 2,
            arg: "x".to_string(),
        };
        assert_eq!(
            derivative_rule(&expr),
            "(2/(x^2 - 1)) * (x*P_2(x) - P_1(x))"
        );
    }

    #[test]
    fn test_derivative_rule_bessel() {
        let expr = PhysicsSpecialExpr::BesselJ {
            n: 1,
            arg: "x".to_string(),
        };
        assert_eq!(derivative_rule(&expr), "(J_0(x) - J_2(x))/2");
    }

    #[test]
    fn test_integral_rule_legendre() {
        let expr = PhysicsSpecialExpr::Legendre {
            n: 2,
            arg: "x".to_string(),
        };
        assert_eq!(
            integral_rule(&expr),
            Some("(P_3(x) - P_1(x))/5 + C".to_string())
        );
    }

    #[test]
    fn test_simplify_rule_special_values() {
        let leg = PhysicsSpecialExpr::Legendre {
            n: 3,
            arg: "1".to_string(),
        };
        let bes = PhysicsSpecialExpr::BesselJ {
            n: 2,
            arg: "0".to_string(),
        };
        assert_eq!(simplify_rule(&leg), Some("1".to_string()));
        assert_eq!(simplify_rule(&bes), Some("0".to_string()));
    }
}
