//! Benchmarks for Phase 2 algorithms.
//!
//! Includes:
//! - Sparse polynomial algorithms (GCD)
//! - Linear algebra (sparse matrix operations)

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use tertius_linalg::sparse_matrix::CsrMatrix;
use tertius_poly::algorithms::sparse_gcd::sparse_gcd_univariate;
use tertius_rings::finite_field::FiniteField;

type GF101 = FiniteField<101>;

fn gf101(n: i64) -> GF101 {
    GF101::new(n.rem_euclid(101) as u64)
}

/// Benchmark sparse GCD.
fn bench_sparse_gcd(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse_gcd");

    // (x+1)(x+2) = x^2 + 3x + 2
    let f = vec![gf101(2), gf101(3), gf101(1)];
    // (x+1)(x+3) = x^2 + 4x + 3
    let g = vec![gf101(3), gf101(4), gf101(1)];

    group.bench_function("degree_2_common_factor", |b| {
        b.iter(|| black_box(sparse_gcd_univariate::<101>(&f, &g)))
    });

    // Larger example
    // f = (x+1)^5 mod 101
    let f_large = vec![gf101(1), gf101(5), gf101(10), gf101(10), gf101(5), gf101(1)];
    // g = (x+1)^3 mod 101
    let g_large = vec![gf101(1), gf101(3), gf101(3), gf101(1)];

    group.bench_function("degree_5_3_common_factor", |b| {
        b.iter(|| black_box(sparse_gcd_univariate::<101>(&f_large, &g_large)))
    });

    group.finish();
}

/// Benchmark sparse matrix operations.
fn bench_sparse_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse_matrix");

    // Create a sparse matrix
    for size in [10, 50, 100] {
        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptrs = vec![0];

        // Create a bidiagonal matrix
        for i in 0..size {
            values.push(gf101(1));
            col_indices.push(i);
            if i + 1 < size {
                values.push(gf101(2));
                col_indices.push(i + 1);
            }
            row_ptrs.push(values.len());
        }

        let matrix = CsrMatrix::from_raw(values, col_indices, row_ptrs, size);
        let vector: Vec<GF101> = (0..size).map(|i| gf101(i as i64)).collect();

        group.bench_with_input(
            BenchmarkId::new("spmv", size),
            &size,
            |b, _| b.iter(|| black_box(matrix.spmv(&vector))),
        );
    }

    group.finish();
}

criterion_group!(
    phase2_benches,
    bench_sparse_gcd,
    bench_sparse_matrix,
);

criterion_main!(phase2_benches);
