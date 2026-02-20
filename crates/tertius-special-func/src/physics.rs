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
}
