//! Separation-of-variables starter templates for core PDEs.

/// A separable mode `X(x) * T(t)` or `X(x) * Y(y)`.
#[derive(Clone, Debug)]
pub struct SeparatedMode {
    /// Text representation for human-readable output.
    pub display: String,
}

/// Fourier-mode coefficient for canonical PDE expansions.
#[derive(Clone, Copy, Debug)]
pub struct ModeCoefficient {
    /// Mode index `n`.
    pub n: u32,
    /// Scalar coefficient for this mode.
    pub amplitude: f64,
}

/// Heat equation mode for `u_t = alpha^2 u_xx` on `[0, L]`.
#[must_use]
pub fn heat_mode(n: u32, length: f64, alpha: f64) -> SeparatedMode {
    let n_f = f64::from(n);
    let lambda = n_f * std::f64::consts::PI / length;
    SeparatedMode {
        display: format!("sin({n}*pi*x/{length}) * exp(-({alpha}^2)*({lambda}^2)*t)"),
    }
}

/// Wave equation mode for `u_tt = c^2 u_xx` on `[0, L]`.
#[must_use]
pub fn wave_mode(n: u32, length: f64, c: f64) -> SeparatedMode {
    let n_f = f64::from(n);
    let omega = c * n_f * std::f64::consts::PI / length;
    SeparatedMode {
        display: format!("sin({n}*pi*x/{length}) * cos({omega}*t)"),
    }
}

/// Laplace equation rectangle mode for `u_xx + u_yy = 0`.
#[must_use]
pub fn laplace_rectangle_mode(n: u32, length_x: f64) -> SeparatedMode {
    SeparatedMode {
        display: format!("sin({n}*pi*x/{length_x}) * sinh({n}*pi*y/{length_x})"),
    }
}

/// Evaluate a finite heat-series solution on `[0, L]`:
/// `u(x,t)=Σ A_n sin(nπx/L) exp(-alpha^2 (nπ/L)^2 t)`.
#[must_use]
pub fn heat_solution_dirichlet(
    coeffs: &[ModeCoefficient],
    length: f64,
    alpha: f64,
    x: f64,
    t: f64,
) -> f64 {
    coeffs
        .iter()
        .filter(|c| c.n > 0)
        .map(|c| {
            let n = f64::from(c.n);
            let lambda = n * std::f64::consts::PI / length;
            c.amplitude * (lambda * x).sin() * (-(alpha * alpha) * lambda * lambda * t).exp()
        })
        .sum()
}

/// Evaluate a finite wave-series solution on `[0, L]`:
/// `u(x,t)=Σ [A_n cos(ω_n t)+B_n sin(ω_n t)] sin(nπx/L)`.
#[must_use]
pub fn wave_solution_dirichlet(
    cos_coeffs: &[ModeCoefficient],
    sin_coeffs: &[ModeCoefficient],
    length: f64,
    c: f64,
    x: f64,
    t: f64,
) -> f64 {
    let cos_sum: f64 = cos_coeffs
        .iter()
        .filter(|coef| coef.n > 0)
        .map(|coef| {
            let n = f64::from(coef.n);
            let lambda = n * std::f64::consts::PI / length;
            let omega = c * lambda;
            coef.amplitude * (omega * t).cos() * (lambda * x).sin()
        })
        .sum();

    let sin_sum: f64 = sin_coeffs
        .iter()
        .filter(|coef| coef.n > 0)
        .map(|coef| {
            let n = f64::from(coef.n);
            let lambda = n * std::f64::consts::PI / length;
            let omega = c * lambda;
            coef.amplitude * (omega * t).sin() * (lambda * x).sin()
        })
        .sum();

    cos_sum + sin_sum
}

/// Evaluate a finite Laplace-rectangle series on `0<y<H`:
/// `u(x,y)=Σ A_n sin(nπx/L) sinh(nπy/L)/sinh(nπH/L)`.
#[must_use]
pub fn laplace_rectangle_solution(
    coeffs: &[ModeCoefficient],
    length_x: f64,
    height_y: f64,
    x: f64,
    y: f64,
) -> f64 {
    coeffs
        .iter()
        .filter(|c| c.n > 0)
        .map(|c| {
            let n = f64::from(c.n);
            let lambda = n * std::f64::consts::PI / length_x;
            let denom = (lambda * height_y).sinh();
            if denom.abs() < 1e-15 {
                0.0
            } else {
                c.amplitude * (lambda * x).sin() * (lambda * y).sinh() / denom
            }
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heat_mode_contains_expected_terms() {
        let mode = heat_mode(1, 2.0, 0.5);
        assert!(mode.display.contains("sin("));
        assert!(mode.display.contains("exp("));
        assert!(mode.display.contains("pi"));
    }

    #[test]
    fn test_wave_mode_contains_expected_terms() {
        let mode = wave_mode(2, 1.0, 3.0);
        assert!(mode.display.contains("sin("));
        assert!(mode.display.contains("cos("));
    }

    #[test]
    fn test_laplace_rectangle_mode_contains_expected_terms() {
        let mode = laplace_rectangle_mode(1, 1.0);
        assert!(mode.display.contains("sin("));
        assert!(mode.display.contains("sinh("));
    }

    #[test]
    fn test_heat_solution_matches_single_mode_initial_shape() {
        let coeffs = [ModeCoefficient {
            n: 1,
            amplitude: 2.0,
        }];
        let u = heat_solution_dirichlet(&coeffs, 2.0, 0.7, 0.5, 0.0);
        let expected = 2.0 * (std::f64::consts::PI * 0.5 / 2.0).sin();
        assert!((u - expected).abs() < 1e-9);
    }

    #[test]
    fn test_wave_solution_uses_cos_and_sin_modes() {
        let cos_coeffs = [ModeCoefficient {
            n: 1,
            amplitude: 1.0,
        }];
        let sin_coeffs = [ModeCoefficient {
            n: 1,
            amplitude: 0.5,
        }];
        let u = wave_solution_dirichlet(&cos_coeffs, &sin_coeffs, 1.0, 2.0, 0.25, 0.3);
        let omega = 2.0 * std::f64::consts::PI;
        let expected =
            ((omega * 0.3).cos() + 0.5 * (omega * 0.3).sin()) * (std::f64::consts::PI * 0.25).sin();
        assert!((u - expected).abs() < 1e-9);
    }

    #[test]
    fn test_laplace_rectangle_solution_satisfies_zero_bottom_boundary() {
        let coeffs = [ModeCoefficient {
            n: 1,
            amplitude: 3.0,
        }];
        let u = laplace_rectangle_solution(&coeffs, 2.0, 1.0, 0.7, 0.0);
        assert!(u.abs() < 1e-12);
    }
}
