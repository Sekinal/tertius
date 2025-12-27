//! Benchmarks for polynomial multiplication algorithms.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use tertius_integers::Integer;
use tertius_poly::algorithms::fft::fft_multiply_integers;
use tertius_poly::algorithms::ntt::{ntt_multiply, NttField};
use tertius_poly::DensePoly;
use tertius_rings::rationals::Q;

/// Generates a random polynomial with integer coefficients.
fn random_poly_q(degree: usize) -> DensePoly<Q> {
    let coeffs: Vec<Q> = (0..=degree)
        .map(|i| Q::from_integer((i as i64 % 100) - 50))
        .collect();
    DensePoly::new(coeffs)
}

/// Generates a random polynomial for NTT.
fn random_poly_ntt(degree: usize) -> Vec<NttField> {
    (0..=degree)
        .map(|i| NttField::new((i as u64) % 1000))
        .collect()
}

/// Generates a random polynomial with Integer coefficients.
fn random_poly_integer(degree: usize) -> Vec<Integer> {
    (0..=degree)
        .map(|i| Integer::new((i as i64 % 100) - 50))
        .collect()
}

fn bench_polynomial_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_mul");

    // Test different sizes
    for size in [16, 64, 256, 1024] {
        let p = random_poly_q(size);
        let q = random_poly_q(size);

        group.bench_with_input(
            BenchmarkId::new("DensePoly<Q>", size),
            &size,
            |b, _| b.iter(|| black_box(p.mul(&q))),
        );
    }

    group.finish();
}

fn bench_ntt(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_mul");

    for size in [16, 64, 256, 1024, 4096] {
        let a = random_poly_ntt(size);
        let b = random_poly_ntt(size);

        group.bench_with_input(
            BenchmarkId::new("NTT", size),
            &size,
            |b_iter, _| b_iter.iter(|| black_box(ntt_multiply(&a, &b))),
        );
    }

    group.finish();
}

fn bench_fft_integers(c: &mut Criterion) {
    let mut group = c.benchmark_group("fft_integer_mul");

    for size in [16, 64, 256, 1024] {
        let a = random_poly_integer(size);
        let b = random_poly_integer(size);

        group.bench_with_input(
            BenchmarkId::new("FFT+CRT", size),
            &size,
            |b_iter, _| b_iter.iter(|| black_box(fft_multiply_integers(&a, &b))),
        );
    }

    group.finish();
}

fn bench_algorithm_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("algorithm_comparison");
    group.sample_size(50);

    // Compare schoolbook vs Karatsuba at crossover point
    let small = random_poly_q(30);
    let medium = random_poly_q(100);
    let large = random_poly_q(500);

    group.bench_function("degree_30", |b| b.iter(|| black_box(small.mul(&small))));
    group.bench_function("degree_100", |b| b.iter(|| black_box(medium.mul(&medium))));
    group.bench_function("degree_500", |b| b.iter(|| black_box(large.mul(&large))));

    group.finish();
}

criterion_group!(
    benches,
    bench_polynomial_multiplication,
    bench_ntt,
    bench_fft_integers,
    bench_algorithm_comparison
);

criterion_main!(benches);
