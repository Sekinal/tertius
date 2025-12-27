//! Integration tests for tertius-linalg.

#[cfg(test)]
mod integration_tests {
    use crate::dense_matrix::DenseMatrix;
    use crate::parallel::{parallel_dot, parallel_rref, ParallelConfig};
    use crate::simd::Gf2x64;
    use crate::sparse_matrix::CsrMatrix;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_sparse_dense_equivalence() {
        // Create a matrix both ways and verify SpMV gives same result
        let dense = vec![
            vec![Q::from_integer(1), Q::from_integer(0), Q::from_integer(2)],
            vec![Q::from_integer(0), Q::from_integer(3), Q::from_integer(0)],
            vec![Q::from_integer(4), Q::from_integer(0), Q::from_integer(5)],
        ];

        let sparse = CsrMatrix::from_dense(&dense);
        let dense_mat = DenseMatrix::from_rows(dense);

        let x = vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)];

        let sparse_result = sparse.spmv(&x);
        let dense_result = dense_mat.mv(&x);

        assert_eq!(sparse_result, dense_result);
    }

    #[test]
    fn test_parallel_vs_sequential() {
        let sparse = CsrMatrix::from_dense(&vec![
            vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)],
            vec![Q::from_integer(4), Q::from_integer(5), Q::from_integer(6)],
            vec![Q::from_integer(7), Q::from_integer(8), Q::from_integer(9)],
        ]);

        let x = vec![Q::from_integer(1), Q::from_integer(1), Q::from_integer(1)];

        let seq_result = sparse.spmv(&x);
        let par_result = sparse.spmv_parallel(&x);

        assert_eq!(seq_result, par_result);
    }

    #[test]
    fn test_gf2_operations() {
        // Test that GF(2) operations work correctly
        let a = Gf2x64::new(0b1010_1010);
        let b = Gf2x64::new(0b1100_1100);

        // Addition is XOR
        let sum = a + b;
        assert_eq!(sum.raw(), 0b0110_0110);

        // Multiplication is AND
        let prod = a * b;
        assert_eq!(prod.raw(), 0b1000_1000);
    }

    #[test]
    fn test_parallel_rref() {
        let mut m = DenseMatrix::from_rows(vec![
            vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(1)],
            vec![Q::from_integer(2), Q::from_integer(4), Q::from_integer(2)],
            vec![Q::from_integer(1), Q::from_integer(3), Q::from_integer(1)],
        ]);

        let config = ParallelConfig {
            parallel_threshold: 1,
            block_size: 2,
        };

        let rank = parallel_rref(&mut m, &config);

        // Matrix has rank 2 (second row is 2x first row)
        assert_eq!(rank, 2);
    }

    #[test]
    fn test_transpose_consistency() {
        let dense = vec![
            vec![Q::from_integer(1), Q::from_integer(2), Q::from_integer(3)],
            vec![Q::from_integer(4), Q::from_integer(5), Q::from_integer(6)],
        ];

        let sparse = CsrMatrix::from_dense(&dense);
        let transposed = sparse.transpose();

        assert_eq!(transposed.num_rows(), 3);
        assert_eq!(transposed.num_cols(), 2);

        // Check specific entries
        assert_eq!(transposed.get(0, 0), Some(&Q::from_integer(1)));
        assert_eq!(transposed.get(0, 1), Some(&Q::from_integer(4)));
        assert_eq!(transposed.get(1, 0), Some(&Q::from_integer(2)));
    }

    #[test]
    fn test_parallel_dot_product() {
        let a = vec![
            Q::from_integer(1),
            Q::from_integer(2),
            Q::from_integer(3),
            Q::from_integer(4),
        ];
        let b = vec![
            Q::from_integer(5),
            Q::from_integer(6),
            Q::from_integer(7),
            Q::from_integer(8),
        ];

        let result = parallel_dot(&a, &b);
        // 1*5 + 2*6 + 3*7 + 4*8 = 5 + 12 + 21 + 32 = 70
        assert_eq!(result, Q::from_integer(70));
    }
}

#[cfg(test)]
mod proptest_tests {
    // Property-based tests will go here
    // Requires proptest strategies for matrices
}
