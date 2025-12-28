//! Benchmarks for power series and limit algorithms.
//!
//! Includes:
//! - Power series arithmetic (add, mul, inverse)
//! - Series composition
//! - Limit computation

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use smallvec::smallvec;

use tertius_core::{ExprArena, ExprNode};
use tertius_core::expr::functions;
use tertius_limits::{compute_limit, Limit};
use tertius_rings::rationals::Q;
use tertius_series::PowerSeries;

fn q(n: i64, d: i64) -> Q {
    Q::new(n, d)
}

/// Benchmark power series addition.
fn bench_series_add(c: &mut Criterion) {
    let mut group = c.benchmark_group("series_add");

    for precision in [10, 50, 100] {
        // Create exp(x) and geometric series
        let exp: PowerSeries<Q> = PowerSeries::exp(precision);
        let geo: PowerSeries<Q> = PowerSeries::geometric(precision);

        group.bench_with_input(
            BenchmarkId::new("exp+geometric", precision),
            &precision,
            |b, _| {
                b.iter(|| {
                    let result = exp.add(&geo);
                    // Force evaluation of first few coefficients
                    for i in 0..5 {
                        black_box(result.coeff(i));
                    }
                })
            },
        );
    }

    group.finish();
}

/// Benchmark power series multiplication (Cauchy product).
fn bench_series_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("series_mul");

    for precision in [10, 20, 50] {
        // (1 + x + x^2/2 + ...) * (1 + x + x^2/2 + ...)
        let exp: PowerSeries<Q> = PowerSeries::exp(precision);

        group.bench_with_input(
            BenchmarkId::new("exp_squared", precision),
            &precision,
            |b, _| {
                b.iter(|| {
                    let result = exp.mul(&exp);
                    // Force evaluation
                    for i in 0..5 {
                        black_box(result.coeff(i));
                    }
                })
            },
        );
    }

    group.finish();
}

/// Benchmark power series inversion.
fn bench_series_inverse(c: &mut Criterion) {
    let mut group = c.benchmark_group("series_inverse");

    for precision in [10, 20, 50] {
        // 1/(1-x) computed via Newton iteration
        let _geo: PowerSeries<Q> = PowerSeries::geometric(precision);

        group.bench_with_input(
            BenchmarkId::new("1/(1-x)", precision),
            &precision,
            |b, _| {
                b.iter(|| {
                    // Invert 1 - x + x^2 - x^3 + ... (alternating signs)
                    let coeffs: Vec<Q> = (0..precision)
                        .map(|i| if i % 2 == 0 { q(1, 1) } else { q(-1, 1) })
                        .collect();
                    let series = PowerSeries::from_coeffs(coeffs);
                    let inv = series.inverse();
                    black_box(inv.map(|s| s.coeff(0)))
                })
            },
        );
    }

    group.finish();
}

/// Benchmark power series derivative/integral.
fn bench_series_calculus(c: &mut Criterion) {
    let mut group = c.benchmark_group("series_calculus");

    let exp: PowerSeries<Q> = PowerSeries::exp(50);

    group.bench_function("derivative_exp", |b| {
        b.iter(|| {
            let deriv = exp.derivative();
            for i in 0..5 {
                black_box(deriv.coeff(i));
            }
        })
    });

    let sin: PowerSeries<Q> = PowerSeries::sin(50);

    group.bench_function("integral_sin", |b| {
        b.iter(|| {
            let integ = sin.integral();
            for i in 0..5 {
                black_box(integ.coeff(i));
            }
        })
    });

    group.finish();
}

/// Benchmark limit computation.
fn bench_limits(c: &mut Criterion) {
    let mut group = c.benchmark_group("limits");

    // lim(x→∞) exp(x)
    group.bench_function("lim_exp_x", |b| {
        b.iter(|| {
            let mut arena = ExprArena::new();
            let x = arena.symbol("x");
            let exp_x = arena.intern(ExprNode::Function {
                id: functions::EXP,
                args: smallvec![x],
            });
            black_box(compute_limit(&mut arena, exp_x, x, Limit::PosInfinity))
        })
    });

    // lim(x→∞) log(x)
    group.bench_function("lim_log_x", |b| {
        b.iter(|| {
            let mut arena = ExprArena::new();
            let x = arena.symbol("x");
            let log_x = arena.intern(ExprNode::Function {
                id: functions::LN,
                args: smallvec![x],
            });
            black_box(compute_limit(&mut arena, log_x, x, Limit::PosInfinity))
        })
    });

    // lim(x→∞) x^2
    group.bench_function("lim_x_squared", |b| {
        b.iter(|| {
            let mut arena = ExprArena::new();
            let x = arena.symbol("x");
            let two = arena.integer(2);
            let x_sq = arena.pow(x, two);
            black_box(compute_limit(&mut arena, x_sq, x, Limit::PosInfinity))
        })
    });

    // lim(x→∞) exp(-x)
    group.bench_function("lim_exp_neg_x", |b| {
        b.iter(|| {
            let mut arena = ExprArena::new();
            let x = arena.symbol("x");
            let neg_one = arena.integer(-1);
            let neg_x = arena.mul(smallvec![neg_one, x]);
            let exp_neg_x = arena.intern(ExprNode::Function {
                id: functions::EXP,
                args: smallvec![neg_x],
            });
            black_box(compute_limit(&mut arena, exp_neg_x, x, Limit::PosInfinity))
        })
    });

    // lim(x→0) x
    group.bench_function("lim_x_to_0", |b| {
        b.iter(|| {
            let mut arena = ExprArena::new();
            let x = arena.symbol("x");
            let zero = arena.integer(0);
            black_box(compute_limit(&mut arena, x, x, Limit::TwoSided(zero)))
        })
    });

    group.finish();
}

/// Benchmark well-known series generation.
fn bench_standard_series(c: &mut Criterion) {
    let mut group = c.benchmark_group("standard_series");

    for precision in [10, 50, 100] {
        group.bench_with_input(
            BenchmarkId::new("exp", precision),
            &precision,
            |b, &p| {
                b.iter(|| {
                    let exp: PowerSeries<Q> = PowerSeries::exp(p);
                    for i in 0..5 {
                        black_box(exp.coeff(i));
                    }
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("sin", precision),
            &precision,
            |b, &p| {
                b.iter(|| {
                    let sin: PowerSeries<Q> = PowerSeries::sin(p);
                    for i in 0..5 {
                        black_box(sin.coeff(i));
                    }
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("cos", precision),
            &precision,
            |b, &p| {
                b.iter(|| {
                    let cos: PowerSeries<Q> = PowerSeries::cos(p);
                    for i in 0..5 {
                        black_box(cos.coeff(i));
                    }
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    series_benches,
    bench_series_add,
    bench_series_mul,
    bench_series_inverse,
    bench_series_calculus,
    bench_limits,
    bench_standard_series,
);

criterion_main!(series_benches);
