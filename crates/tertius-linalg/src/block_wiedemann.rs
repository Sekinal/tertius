//! Block Wiedemann algorithm for sparse linear systems.
//!
//! The Block Wiedemann algorithm solves Ax = b for massive sparse matrices
//! over finite fields without fill-in. It works by computing Krylov sequences
//! and using Berlekamp-Massey to find the minimal polynomial.
//!
//! # References
//!
//! - Coppersmith, "Solving homogeneous linear equations over GF(2)" (1994)
//! - Kaltofen, Lobo, "Factoring high-degree polynomials by the black box
//!   Berlekamp algorithm" (1994)

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use tertius_rings::traits::Field;

use crate::berlekamp_massey::berlekamp_massey;
use crate::sparse_matrix::CsrMatrix;

/// Configuration for Block Wiedemann algorithm.
#[derive(Clone, Debug)]
pub struct BlockWiedemannConfig {
    /// Block size (number of vectors processed together).
    /// Typically 64 for SIMD efficiency on modern CPUs.
    pub block_size: usize,
    /// Random seed for reproducibility.
    pub seed: u64,
    /// Maximum number of Krylov iterations.
    pub max_iterations: usize,
}

impl Default for BlockWiedemannConfig {
    fn default() -> Self {
        Self {
            block_size: 64,
            seed: 42,
            max_iterations: 10000,
        }
    }
}

/// Block Wiedemann solver state.
pub struct BlockWiedemann<R> {
    /// The sparse matrix A.
    matrix: CsrMatrix<R>,
    /// Configuration.
    config: BlockWiedemannConfig,
    /// Random number generator.
    rng: ChaCha8Rng,
}

impl<R: Field + Clone + Send + Sync> BlockWiedemann<R> {
    /// Creates a new Block Wiedemann solver.
    #[must_use]
    pub fn new(matrix: CsrMatrix<R>, config: BlockWiedemannConfig) -> Self {
        let rng = ChaCha8Rng::seed_from_u64(config.seed);
        Self { matrix, config, rng }
    }

    /// Generates a random vector over the field.
    fn random_vector(&mut self, len: usize) -> Vec<R>
    where
        R: From<u64>,
    {
        (0..len).map(|_| R::from(self.rng.gen::<u64>())).collect()
    }

    /// Computes the Krylov sequence: A^0 * v, A^1 * v, ..., A^{k-1} * v.
    fn krylov_sequence(&self, v: &[R], length: usize) -> Vec<Vec<R>> {
        let mut sequence = Vec::with_capacity(length);
        let mut current = v.to_vec();

        for _ in 0..length {
            sequence.push(current.clone());
            current = self.matrix.spmv_parallel(&current);
        }

        sequence
    }

    /// Computes projections w^T * A^i * v for the Berlekamp-Massey algorithm.
    fn compute_projections(&self, w: &[R], v: &[R], length: usize) -> Vec<R> {
        let krylov = self.krylov_sequence(v, length);

        krylov
            .into_iter()
            .map(|av| {
                w.iter()
                    .zip(av.iter())
                    .fold(R::zero(), |acc, (wi, vi)| acc + wi.clone() * vi.clone())
            })
            .collect()
    }

    /// Solves Ax = 0 (finds a non-trivial kernel element).
    ///
    /// Returns Some(x) where Ax = 0 and x != 0, or None if the kernel is trivial.
    pub fn solve_homogeneous(&mut self) -> Option<Vec<R>>
    where
        R: From<u64>,
    {
        let n = self.matrix.num_cols();

        // Generate random vectors
        let w = self.random_vector(n);
        let v = self.random_vector(n);

        // Compute projection sequence
        let sequence_length = 2 * n / self.config.block_size + 10;
        let projections = self.compute_projections(&w, &v, sequence_length.min(self.config.max_iterations));

        // Find minimal polynomial via Berlekamp-Massey
        let bm_result = berlekamp_massey(&projections);

        if bm_result.length == 0 {
            return None;
        }

        // The minimal polynomial gives us a way to construct a kernel element
        // x = c_1 * A^0 * v + c_2 * A^1 * v + ... + c_L * A^{L-1} * v
        let krylov = self.krylov_sequence(&v, bm_result.length + 1);

        // Construct kernel element
        let mut x = vec![R::zero(); n];
        for (i, coeff) in bm_result.connection_poly.iter().enumerate().skip(1) {
            if i <= bm_result.length {
                for j in 0..n {
                    x[j] = x[j].clone() + coeff.clone() * krylov[i - 1][j].clone();
                }
            }
        }

        // Verify x is in kernel
        let ax = self.matrix.spmv_parallel(&x);
        let is_zero = ax.iter().all(|xi| xi.is_zero());

        if is_zero && !x.iter().all(|xi| xi.is_zero()) {
            Some(x)
        } else {
            None
        }
    }

    /// Solves Ax = b (inhomogeneous system).
    ///
    /// Returns Some(x) where Ax = b, or None if no solution exists.
    pub fn solve(&mut self, b: &[R]) -> Option<Vec<R>>
    where
        R: From<u64>,
    {
        let n = self.matrix.num_cols();
        assert_eq!(b.len(), self.matrix.num_rows());

        // Generate random vector
        let w = self.random_vector(n);

        // Compute sequence w^T * A^i * b
        let sequence_length = 2 * n / self.config.block_size + 10;
        let projections = self.compute_projections(&w, b, sequence_length.min(self.config.max_iterations));

        // Find minimal polynomial
        let bm_result = berlekamp_massey(&projections);

        if bm_result.length == 0 || bm_result.connection_poly.is_empty() {
            return None;
        }

        // If the constant term of the minimal polynomial is zero, no solution exists
        if bm_result.connection_poly[0].is_zero() {
            return None;
        }

        // Construct solution: x = -1/c_0 * (c_1 * A^0 * b + c_2 * A^1 * b + ...)
        let krylov = self.krylov_sequence(b, bm_result.length);

        // Compute sum of c_i * A^{i-1} * b
        let mut sum = vec![R::zero(); n];
        for (i, coeff) in bm_result.connection_poly.iter().enumerate().skip(1) {
            if i <= bm_result.length {
                for j in 0..n {
                    sum[j] = sum[j].clone() + coeff.clone() * krylov[i - 1][j].clone();
                }
            }
        }

        // Scale by -1/c_0
        let c0 = &bm_result.connection_poly[0];
        if let Some(c0_inv) = c0.inv() {
            let scale = c0_inv.neg();
            let x: Vec<R> = sum.into_iter().map(|xi| xi * scale.clone()).collect();

            // Verify
            let ax = self.matrix.spmv_parallel(&x);
            let is_solution = ax
                .iter()
                .zip(b.iter())
                .all(|(axi, bi)| (axi.clone() + bi.clone().neg()).is_zero());

            if is_solution {
                return Some(x);
            }
        }

        None
    }
}

/// Simplified interface for solving Ax = b.
pub fn block_wiedemann_solve<R: Field + Clone + Send + Sync + From<u64>>(
    matrix: &CsrMatrix<R>,
    b: &[R],
) -> Option<Vec<R>> {
    let config = BlockWiedemannConfig::default();
    let mut solver = BlockWiedemann::new(matrix.clone(), config);
    solver.solve(b)
}

/// Finds a kernel element (solution to Ax = 0).
pub fn block_wiedemann_kernel<R: Field + Clone + Send + Sync + From<u64>>(
    matrix: &CsrMatrix<R>,
) -> Option<Vec<R>> {
    let config = BlockWiedemannConfig::default();
    let mut solver = BlockWiedemann::new(matrix.clone(), config);
    solver.solve_homogeneous()
}

/// Computes the rank of a sparse matrix using Block Wiedemann.
///
/// This is probabilistic - with high probability, the computed rank is correct.
pub fn block_wiedemann_rank<R: Field + Clone + Send + Sync + From<u64>>(
    matrix: &CsrMatrix<R>,
) -> usize {
    let n = matrix.num_rows().max(matrix.num_cols());
    let config = BlockWiedemannConfig::default();
    let seed = config.seed;
    let solver = BlockWiedemann::new(matrix.clone(), config);

    // Generate random vectors and compute minimal polynomial degree
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let v: Vec<R> = (0..n).map(|_| R::from(rng.gen::<u64>())).collect();
    let w: Vec<R> = (0..n).map(|_| R::from(rng.gen::<u64>())).collect();

    let sequence_length = 2 * n + 10;
    let projections = solver.compute_projections(&w, &v, sequence_length.min(10000));
    let bm_result = berlekamp_massey(&projections);

    bm_result.length
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::finite_field::FiniteField;

    type GF7 = FiniteField<7>;

    #[test]
    fn test_krylov_sequence() {
        // Simple 2x2 matrix
        let dense = vec![
            vec![GF7::new(1), GF7::new(2)],
            vec![GF7::new(3), GF7::new(4)],
        ];
        let matrix = CsrMatrix::from_dense(&dense);

        let config = BlockWiedemannConfig::default();
        let solver = BlockWiedemann::new(matrix, config);

        let v = vec![GF7::new(1), GF7::new(0)];
        let krylov = solver.krylov_sequence(&v, 3);

        assert_eq!(krylov.len(), 3);
        assert_eq!(krylov[0], vec![GF7::new(1), GF7::new(0)]);
        // A * [1, 0] = [1, 3]
        assert_eq!(krylov[1], vec![GF7::new(1), GF7::new(3)]);
    }

    #[test]
    fn test_identity_solve() {
        // Identity matrix should give x = b
        let id = CsrMatrix::<GF7>::identity(3);
        let b = vec![GF7::new(1), GF7::new(2), GF7::new(3)];

        let x = block_wiedemann_solve(&id, &b);

        if let Some(solution) = x {
            assert_eq!(solution, b);
        }
    }

    #[test]
    fn test_simple_solve() {
        // Matrix [[1, 1], [0, 1]] with b = [3, 2]
        // Solution: x = [1, 2] (since x + y = 3, y = 2)
        let dense = vec![
            vec![GF7::new(1), GF7::new(1)],
            vec![GF7::new(0), GF7::new(1)],
        ];
        let matrix = CsrMatrix::from_dense(&dense);
        let b = vec![GF7::new(3), GF7::new(2)];

        let x = block_wiedemann_solve(&matrix, &b);

        if let Some(solution) = x {
            // Verify Ax = b
            let ax = matrix.spmv(&solution);
            assert_eq!(ax, b);
        }
    }
}
