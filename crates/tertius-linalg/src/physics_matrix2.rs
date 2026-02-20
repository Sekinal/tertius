//! Small 2x2 real matrix utilities for physics workflows.

/// Real 2x2 matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix2 {
    /// Entry (1,1)
    pub a11: f64,
    /// Entry (1,2)
    pub a12: f64,
    /// Entry (2,1)
    pub a21: f64,
    /// Entry (2,2)
    pub a22: f64,
}

/// Real eigen-decomposition data for a 2x2 matrix with distinct real eigenvalues.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealEigenDecomposition2 {
    /// First eigenvalue.
    pub lambda1: f64,
    /// First normalized eigenvector.
    pub v1: [f64; 2],
    /// Second eigenvalue.
    pub lambda2: f64,
    /// Second normalized eigenvector.
    pub v2: [f64; 2],
}

impl Matrix2 {
    /// Construct a matrix.
    #[must_use]
    pub fn new(a11: f64, a12: f64, a21: f64, a22: f64) -> Self {
        Self { a11, a12, a21, a22 }
    }

    /// Trace.
    #[must_use]
    pub fn trace(self) -> f64 {
        self.a11 + self.a22
    }

    /// Determinant.
    #[must_use]
    pub fn det(self) -> f64 {
        self.a11 * self.a22 - self.a12 * self.a21
    }

    /// Real eigenvalues if discriminant is non-negative.
    #[must_use]
    pub fn eigenvalues_real(self) -> Option<(f64, f64)> {
        let tr = self.trace();
        let disc = tr * tr - 4.0 * self.det();
        if disc < 0.0 {
            return None;
        }
        let s = disc.sqrt();
        Some(((tr - s) / 2.0, (tr + s) / 2.0))
    }

    /// Matrix exponential exp(A).
    #[must_use]
    pub fn exp(self) -> Self {
        let mu = self.trace() / 2.0;
        let e_mu = mu.exp();

        let b11 = self.a11 - mu;
        let b12 = self.a12;
        let b21 = self.a21;
        let b22 = self.a22 - mu;

        // For a 2x2 trace-zero matrix B, B^2 = delta2 * I.
        let delta2 = b11 * b11 + b12 * b21;
        let (c0, c1) = if delta2 > 1e-15 {
            let d = delta2.sqrt();
            (d.cosh(), d.sinh() / d)
        } else if delta2 < -1e-15 {
            let w = (-delta2).sqrt();
            (w.cos(), w.sin() / w)
        } else {
            (1.0, 1.0)
        };

        let m11 = e_mu * (c0 + c1 * b11);
        let m12 = e_mu * (c1 * b12);
        let m21 = e_mu * (c1 * b21);
        let m22 = e_mu * (c0 + c1 * b22);
        Self::new(m11, m12, m21, m22)
    }

    /// Applies matrix to a 2D vector.
    #[must_use]
    pub fn apply(self, v: [f64; 2]) -> [f64; 2] {
        [
            self.a11 * v[0] + self.a12 * v[1],
            self.a21 * v[0] + self.a22 * v[1],
        ]
    }

    /// Returns `exp(self * t)`.
    #[must_use]
    pub fn state_transition(self, t: f64) -> Self {
        Self::new(self.a11 * t, self.a12 * t, self.a21 * t, self.a22 * t).exp()
    }

    /// Real modal decomposition for distinct real eigenvalues.
    #[must_use]
    pub fn real_modal_decomposition(self) -> Option<RealEigenDecomposition2> {
        let (l1, l2) = self.eigenvalues_real()?;
        if (l1 - l2).abs() < 1e-12 {
            return None;
        }
        let v1 = eigenvector_for(self, l1);
        let v2 = eigenvector_for(self, l2);
        Some(RealEigenDecomposition2 {
            lambda1: l1,
            v1,
            lambda2: l2,
            v2,
        })
    }
}

fn eigenvector_for(m: Matrix2, lambda: f64) -> [f64; 2] {
    let m11 = m.a11 - lambda;
    let m22 = m.a22 - lambda;

    let candidate_a = [m.a12, -m11];
    let candidate_b = [-m22, m.a21];

    let norm_a = (candidate_a[0] * candidate_a[0] + candidate_a[1] * candidate_a[1]).sqrt();
    let norm_b = (candidate_b[0] * candidate_b[0] + candidate_b[1] * candidate_b[1]).sqrt();

    let v = if norm_a >= norm_b {
        candidate_a
    } else {
        candidate_b
    };
    let norm = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if norm < 1e-12 {
        [1.0, 0.0]
    } else {
        [v[0] / norm, v[1] / norm]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eigenvalues_diagonal() {
        let m = Matrix2::new(2.0, 0.0, 0.0, 3.0);
        let (l1, l2) = m.eigenvalues_real().unwrap();
        assert!((l1 - 2.0).abs() < 1e-12 || (l1 - 3.0).abs() < 1e-12);
        assert!((l2 - 2.0).abs() < 1e-12 || (l2 - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_matrix_exponential_rotation_generator() {
        // A = [[0, 1], [-1, 0]] => exp(A) = [[cos(1), sin(1)], [-sin(1), cos(1)]]
        let a = Matrix2::new(0.0, 1.0, -1.0, 0.0);
        let e = a.exp();
        let c = 1.0f64.cos();
        let s = 1.0f64.sin();
        assert!((e.a11 - c).abs() < 1e-9);
        assert!((e.a12 - s).abs() < 1e-9);
        assert!((e.a21 + s).abs() < 1e-9);
        assert!((e.a22 - c).abs() < 1e-9);
    }

    #[test]
    fn test_real_modal_decomposition_diagonal() {
        let m = Matrix2::new(2.0, 0.0, 0.0, 5.0);
        let modal = m.real_modal_decomposition().unwrap();
        assert!((modal.lambda1 - 2.0).abs() < 1e-12 || (modal.lambda2 - 2.0).abs() < 1e-12);
        assert!((modal.lambda1 - 5.0).abs() < 1e-12 || (modal.lambda2 - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_state_transition_diagonal_system() {
        let a = Matrix2::new(-1.0, 0.0, 0.0, -2.0);
        let phi = a.state_transition(1.0);
        let y = phi.apply([1.0, 2.0]);
        assert!((y[0] - (-1.0f64).exp()).abs() < 1e-9);
        assert!((y[1] - 2.0 * (-2.0f64).exp()).abs() < 1e-9);
    }
}
