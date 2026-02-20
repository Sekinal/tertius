//! Separation-of-variables starter templates for core PDEs.

/// A separable mode `X(x) * T(t)` or `X(x) * Y(y)`.
#[derive(Clone, Debug)]
pub struct SeparatedMode {
    /// Text representation for human-readable output.
    pub display: String,
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
}
