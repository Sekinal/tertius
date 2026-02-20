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
}
