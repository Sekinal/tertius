//! Benchmarks for Phase 2 algorithms.
//!
//! Includes:
//! - Sparse polynomial algorithms (GCD)
//! - Linear algebra (sparse matrix operations)
//! - Polynomial factorization (Van Hoeij)
//! - Gröbner basis (M5GB)

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use tertius_factor::van_hoeij_factor;
use tertius_groebner::monomial::PackedMonomial;
use tertius_groebner::M5GB;
use tertius_groebner::m5gb::M5GBConfig;
use tertius_integers::Integer;
use tertius_linalg::sparse_matrix::CsrMatrix;
use tertius_poly::algorithms::sparse_gcd::sparse_gcd_univariate;
use tertius_poly::dense::DensePoly;
use tertius_rings::finite_field::FiniteField;
use tertius_rings::integers::Z;

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

fn z(n: i64) -> Z {
    Z(Integer::new(n))
}

fn poly_z(coeffs: &[i64]) -> DensePoly<Z> {
    DensePoly::new(coeffs.iter().map(|&n| z(n)).collect())
}

/// Benchmark Van Hoeij factorization.
fn bench_factorization(c: &mut Criterion) {
    let mut group = c.benchmark_group("factorization");

    // x^4 - 1 = (x-1)(x+1)(x^2+1)
    let x4_minus_1 = poly_z(&[-1, 0, 0, 0, 1]);
    group.bench_function("x^4-1", |b| {
        b.iter(|| black_box(van_hoeij_factor(&x4_minus_1)))
    });

    // x^6 - 1
    let x6_minus_1 = poly_z(&[-1, 0, 0, 0, 0, 0, 1]);
    group.bench_function("x^6-1", |b| {
        b.iter(|| black_box(van_hoeij_factor(&x6_minus_1)))
    });

    // (x+1)^5
    let x_plus_1_pow_5 = poly_z(&[1, 5, 10, 10, 5, 1]);
    group.bench_function("(x+1)^5", |b| {
        b.iter(|| black_box(van_hoeij_factor(&x_plus_1_pow_5)))
    });

    // x^4 + 4 (Sophie Germain identity)
    let sophie = poly_z(&[4, 0, 0, 0, 1]);
    group.bench_function("x^4+4", |b| {
        b.iter(|| black_box(van_hoeij_factor(&sophie)))
    });

    group.finish();
}

fn ff(n: u64) -> GF101 {
    GF101::new(n % 101)
}

fn mono(exps: &[u16]) -> PackedMonomial {
    PackedMonomial::new(exps)
}

/// Benchmark M5GB Gröbner basis.
fn bench_groebner(c: &mut Criterion) {
    let mut group = c.benchmark_group("groebner");

    // cyclic-3: x + y + z, xy + yz + zx, xyz - 1
    let cyclic3 = vec![
        // x + y + z
        vec![(ff(1), mono(&[1, 0, 0])), (ff(1), mono(&[0, 1, 0])), (ff(1), mono(&[0, 0, 1]))],
        // xy + yz + zx
        vec![(ff(1), mono(&[1, 1, 0])), (ff(1), mono(&[0, 1, 1])), (ff(1), mono(&[1, 0, 1]))],
        // xyz - 1
        vec![(ff(1), mono(&[1, 1, 1])), (ff(100), mono(&[0, 0, 0]))],  // 100 = -1 mod 101
    ];

    group.bench_function("cyclic-3", |b| {
        b.iter(|| {
            let m5gb = M5GB::new(cyclic3.clone(), M5GBConfig::default());
            black_box(m5gb.compute())
        })
    });

    // katsura-3: x + 2y + 2z - 1, x^2 + 2y^2 + 2z^2 - x, 2xy + 2yz - y
    let katsura3 = vec![
        // x + 2y + 2z - 1
        vec![
            (ff(1), mono(&[1, 0, 0])),
            (ff(2), mono(&[0, 1, 0])),
            (ff(2), mono(&[0, 0, 1])),
            (ff(100), mono(&[0, 0, 0])),
        ],
        // x^2 + 2y^2 + 2z^2 - x
        vec![
            (ff(1), mono(&[2, 0, 0])),
            (ff(2), mono(&[0, 2, 0])),
            (ff(2), mono(&[0, 0, 2])),
            (ff(100), mono(&[1, 0, 0])),
        ],
        // 2xy + 2yz - y
        vec![
            (ff(2), mono(&[1, 1, 0])),
            (ff(2), mono(&[0, 1, 1])),
            (ff(100), mono(&[0, 1, 0])),
        ],
    ];

    group.bench_function("katsura-3", |b| {
        b.iter(|| {
            let m5gb = M5GB::new(katsura3.clone(), M5GBConfig::default());
            black_box(m5gb.compute())
        })
    });

    group.finish();
}

criterion_group!(
    phase2_benches,
    bench_sparse_gcd,
    bench_sparse_matrix,
    bench_factorization,
    bench_groebner,
);

criterion_main!(phase2_benches);
