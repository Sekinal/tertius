//! ODE starter solvers for physics workflows.

/// Closed-form first-order linear ODE solution.
#[derive(Clone, Debug)]
pub struct FirstOrderLinearSolution {
    a_coeff: f64,
    b: f64,
    c_integration: f64,
}

impl FirstOrderLinearSolution {
    /// Evaluate y(x).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        if self.a_coeff.abs() < 1e-12 {
            self.b * x + self.c_integration
        } else {
            let steady = self.b / self.a_coeff;
            steady + self.c_integration * (-self.a_coeff * x).exp()
        }
    }
}

/// Solve `y' + a y = b` with optional initial value `y(0)=y0`.
#[must_use]
pub fn solve_first_order_linear_constant(
    a: f64,
    b: f64,
    y0: Option<f64>,
) -> FirstOrderLinearSolution {
    let c_integration = if a.abs() < 1e-12 {
        y0.unwrap_or(0.0)
    } else {
        y0.map_or(0.0, |y_init| y_init - b / a)
    };
    FirstOrderLinearSolution {
        a_coeff: a,
        b,
        c_integration,
    }
}

/// Closed-form second-order homogeneous constant-coefficient solution.
#[derive(Clone, Debug)]
pub struct SecondOrderHomogeneousSolution {
    alpha: f64,
    beta: f64,
    c1: f64,
    c2: f64,
    mode: SecondOrderMode,
}

#[derive(Clone, Copy, Debug)]
enum SecondOrderMode {
    DistinctReal { r1: f64, r2: f64 },
    RepeatedReal { r: f64 },
    ComplexPair,
}

impl SecondOrderHomogeneousSolution {
    /// Evaluate y(x).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        match self.mode {
            SecondOrderMode::DistinctReal { r1, r2 } => {
                self.c1 * (r1 * x).exp() + self.c2 * (r2 * x).exp()
            }
            SecondOrderMode::RepeatedReal { r } => (self.c1 + self.c2 * x) * (r * x).exp(),
            SecondOrderMode::ComplexPair => {
                let envelope = (self.alpha * x).exp();
                envelope * (self.c1 * (self.beta * x).cos() + self.c2 * (self.beta * x).sin())
            }
        }
    }
}

/// Solve `a y'' + b y' + c y = 0` with optional initial conditions `(y(0), y'(0))`.
#[must_use]
pub fn solve_second_order_constant_homogeneous(
    a: f64,
    b: f64,
    c: f64,
    initial: Option<(f64, f64)>,
) -> SecondOrderHomogeneousSolution {
    if a.abs() < 1e-12 {
        // Degenerate fallback: by' + cy = 0
        let lambda = if b.abs() < 1e-12 { 0.0 } else { -c / b };
        let (y0, _) = initial.unwrap_or((0.0, 0.0));
        return SecondOrderHomogeneousSolution {
            alpha: lambda,
            beta: 0.0,
            c1: y0,
            c2: 0.0,
            mode: SecondOrderMode::DistinctReal {
                r1: lambda,
                r2: lambda,
            },
        };
    }

    let (y0, v0) = initial.unwrap_or((0.0, 0.0));
    let disc = b * b - 4.0 * a * c;
    if disc > 1e-12 {
        let root_disc = disc.sqrt();
        let r1 = (-b + root_disc) / (2.0 * a);
        let r2 = (-b - root_disc) / (2.0 * a);
        let c1 = (v0 - y0 * r2) / (r1 - r2);
        let c2 = y0 - c1;
        return SecondOrderHomogeneousSolution {
            alpha: 0.0,
            beta: 0.0,
            c1,
            c2,
            mode: SecondOrderMode::DistinctReal { r1, r2 },
        };
    }

    if disc.abs() <= 1e-12 {
        let r = -b / (2.0 * a);
        let c1 = y0;
        let c2 = v0 - r * y0;
        return SecondOrderHomogeneousSolution {
            alpha: 0.0,
            beta: 0.0,
            c1,
            c2,
            mode: SecondOrderMode::RepeatedReal { r },
        };
    }

    // Complex conjugate roots alpha ± i beta.
    let alpha = -b / (2.0 * a);
    let beta = (-disc).sqrt() / (2.0 * a);
    let c1 = y0;
    let c2 = (v0 - alpha * y0) / beta;
    SecondOrderHomogeneousSolution {
        alpha,
        beta,
        c1,
        c2,
        mode: SecondOrderMode::ComplexPair,
    }
}

/// Closed-form separable ODE solution for `y' = k y^n`.
#[derive(Clone, Debug)]
pub struct SeparablePowerSolution {
    k: f64,
    n: f64,
    x0: f64,
    y0: f64,
}

impl SeparablePowerSolution {
    /// Evaluate y(x).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        if (self.n - 1.0).abs() < 1e-12 {
            self.y0 * (self.k * (x - self.x0)).exp()
        } else {
            let one_minus_n = 1.0 - self.n;
            let base = self.y0.powf(one_minus_n) + one_minus_n * self.k * (x - self.x0);
            base.powf(1.0 / one_minus_n)
        }
    }
}

/// Solve `y' = k y^n` with initial condition `y(x0)=y0`.
#[must_use]
pub fn solve_separable_power(k: f64, n: f64, x0: f64, y0: f64) -> SeparablePowerSolution {
    SeparablePowerSolution { k, n, x0, y0 }
}

/// Residual `y' + a y - b` using centered finite differences.
#[must_use]
pub fn verify_first_order_linear_constant_solution(
    solution: &FirstOrderLinearSolution,
    a: f64,
    b: f64,
    x: f64,
) -> f64 {
    let h = 1e-6;
    let y_plus = solution.eval(x + h);
    let y_minus = solution.eval(x - h);
    let y = solution.eval(x);
    let dy = (y_plus - y_minus) / (2.0 * h);
    dy + a * y - b
}

/// Residual `a y'' + b y' + c y` using finite differences.
#[must_use]
pub fn verify_second_order_constant_homogeneous_solution(
    solution: &SecondOrderHomogeneousSolution,
    a: f64,
    b: f64,
    c: f64,
    x: f64,
) -> f64 {
    let h = 1e-5;
    let y_plus = solution.eval(x + h);
    let y = solution.eval(x);
    let y_minus = solution.eval(x - h);
    let dy = (y_plus - y_minus) / (2.0 * h);
    let ddy = (y_plus - 2.0 * y + y_minus) / (h * h);
    a * ddy + b * dy + c * y
}

/// Numerical fallback by forward Euler for `y' + a y = b`.
#[must_use]
pub fn euler_fallback_first_order_linear(
    a: f64,
    b: f64,
    y0: f64,
    x0: f64,
    x_target: f64,
    step: f64,
) -> f64 {
    if step.abs() < 1e-15 || (x_target - x0).abs() < 1e-15 {
        return y0;
    }

    let mut x = x0;
    let mut y = y0;
    let h = if x_target >= x0 {
        step.abs()
    } else {
        -step.abs()
    };

    while (x_target - x).abs() > h.abs() {
        y += h * (b - a * y);
        x += h;
    }

    let rem = x_target - x;
    if rem.abs() > 1e-15 {
        y += rem * (b - a * y);
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_first_order_linear_constant_ivp() {
        // y' + 2y = 4, y(0)=1 => y = 2 - e^{-2x}
        let sol = solve_first_order_linear_constant(2.0, 4.0, Some(1.0));
        assert!((sol.eval(0.0) - 1.0).abs() < 1e-9);
        let expected = 2.0 - (-2.0f64).exp();
        assert!((sol.eval(1.0) - expected).abs() < 1e-9);
    }

    #[test]
    fn test_second_order_constant_homogeneous_ivp() {
        // y'' + y = 0, y(0)=0, y'(0)=1 => y = sin(x)
        let sol = solve_second_order_constant_homogeneous(1.0, 0.0, 1.0, Some((0.0, 1.0)));
        assert!(sol.eval(0.0).abs() < 1e-9);
        assert!((sol.eval(1.0) - 1.0f64.sin()).abs() < 1e-9);
    }

    #[test]
    fn test_separable_power_solution() {
        // y' = 2 y^2, y(0)=1 => y = 1/(1 - 2x)
        let sol = solve_separable_power(2.0, 2.0, 0.0, 1.0);
        assert!((sol.eval(0.0) - 1.0).abs() < 1e-9);
        assert!((sol.eval(0.1) - 1.25).abs() < 1e-9);
    }

    #[test]
    fn test_verify_first_order_solution_residual() {
        let sol = solve_first_order_linear_constant(2.0, 4.0, Some(1.0));
        let residual = verify_first_order_linear_constant_solution(&sol, 2.0, 4.0, 0.7);
        assert!(residual.abs() < 1e-5);
    }

    #[test]
    fn test_verify_second_order_solution_residual() {
        let sol = solve_second_order_constant_homogeneous(1.0, 0.0, 1.0, Some((0.0, 1.0)));
        let residual = verify_second_order_constant_homogeneous_solution(&sol, 1.0, 0.0, 1.0, 0.8);
        assert!(residual.abs() < 1e-4);
    }

    #[test]
    fn test_euler_fallback_first_order_linear() {
        // y' + y = 1, y(0)=0 => y(1) = 1 - e^-1
        let y_num = euler_fallback_first_order_linear(1.0, 1.0, 0.0, 0.0, 1.0, 1e-3);
        let y_exact = 1.0 - (-1.0f64).exp();
        assert!((y_num - y_exact).abs() < 2e-3);
    }
}
