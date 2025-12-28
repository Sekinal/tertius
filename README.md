# Tertius

**A next-generation Computer Algebra System in Rust**

Tertius is a high-performance, modular CAS designed for speed and correctness. It leverages Rust's type system and modern algorithms to achieve performance rivaling commercial systems.

[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.75+-orange.svg)](https://www.rust-lang.org/)

## Features

### Core Capabilities (Phase 1)

- **Arbitrary Precision Arithmetic**: Integers and rationals with no size limits
- **Polynomial Algebra**: Dense and sparse representations with adaptive algorithms
- **Modular Arithmetic**: Compile-time and runtime prime fields
- **Algebraic Number Fields**: Q(α) for algebraic extensions
- **Symbolic Simplification**: Equality saturation via the `egg` library

### Advanced Algorithms (Phase 2)

- **Polynomial Factorization**: Van Hoeij algorithm with LLL lattice reduction
- **Gröbner Bases**: M5GB algorithm (F5 signatures + M4GB caching)
- **Sparse Linear Algebra**: Block Wiedemann, Smith Normal Form
- **Sparse Interpolation**: Ben-Or/Tiwari algorithm
- **Sparse GCD**: Hu-Monagan parallel algorithm
- **Symbolic Integration**: Differentiation/integration rules via equality saturation

### Performance Highlights

- **Hash-Consing**: O(1) structural equality via arena allocation
- **Algorithm Selection**: Automatic switching between schoolbook, Karatsuba, and FFT/NTT
- **SIMD-Friendly**: Bit-packed monomials for cache-efficient operations
- **Parallel Throughout**: `rayon`-based parallelism in all major algorithms

## Benchmark Results

### NTT Performance (Modular Arithmetic)

| Degree | Tertius NTT | SymPy ZZ | Speedup |
|--------|-------------|----------|---------|
| 16 | 7 µs | 39 µs | **5.6x** |
| 64 | 25 µs | 255 µs | **10.2x** |
| 256 | 107 µs | 2.8 ms | **26x** |
| 1024 | 512 µs | 22.9 ms | **45x** |
| 4096 | 2.35 ms | - | - |

### Algorithm Comparison

| Operation | Degree | Time | Algorithm |
|-----------|--------|------|-----------|
| Dense Poly (Q) | 30 | 163 µs | Schoolbook |
| Dense Poly (Q) | 100 | 1.36 ms | Karatsuba |
| Dense Poly (Q) | 500 | 14.7 ms | Karatsuba |
| FFT+CRT (Integer) | 256 | 2.55 ms | 3-prime NTT |
| FFT+CRT (Integer) | 1024 | 10.4 ms | 3-prime NTT |

### Phase 2 Algorithms

| Operation | Input Size | Time | Notes |
|-----------|------------|------|-------|
| Sparse GCD | degree 2 | 102 ns | Common factor detection |
| Sparse GCD | degree 5 | 96 ns | Euclidean algorithm over GF(101) |
| Sparse SpMV | 10×10 | 57 ns | Bidiagonal CSR matrix |
| Sparse SpMV | 50×50 | 272 ns | Bidiagonal CSR matrix |
| Sparse SpMV | 100×100 | 547 ns | Bidiagonal CSR matrix |

## Architecture

```
tertius/
├── crates/
│   ├── tertius-core/      # Arena, expressions, hash-consing
│   ├── tertius-integers/  # Arbitrary precision (dashu wrapper)
│   ├── tertius-rings/     # Ring/Field traits, Z, Q, Z_p
│   ├── tertius-poly/      # Polynomial arithmetic + sparse algorithms
│   ├── tertius-simplify/  # egg-based simplification + calculus
│   ├── tertius-linalg/    # Sparse linear algebra (CSR, Wiedemann)
│   ├── tertius-factor/    # Polynomial factorization (Van Hoeij, LLL)
│   ├── tertius-groebner/  # Gröbner bases (M5GB algorithm)
│   └── tertius/           # Facade crate
├── benches/               # Criterion benchmarks
└── scripts/               # Comparison scripts
```

## Quick Start

```rust
use tertius::prelude::*;

// Polynomial arithmetic
let x = DensePoly::<Q>::x();
let one = DensePoly::one();

let p = x.mul(&x).add(&x).add(&one);  // x² + x + 1
let q = x.sub(&one);                   // x - 1

let product = p.mul(&q);
println!("({}) * ({}) = {}", p, q, product);

// Modular arithmetic
use tertius_integers::ModInt;

const P: u64 = 998244353;  // NTT-friendly prime
let a = ModInt::<P>::new(42);
let b = a.inv().unwrap();
assert_eq!((a * b).value(), 1);
```

## Algorithms

### Polynomial Multiplication

| Algorithm | Complexity | When Used |
|-----------|------------|-----------|
| Schoolbook | O(n²) | degree < 32 |
| Karatsuba | O(n^1.58) | 32 ≤ degree < 1024 |
| NTT/FFT | O(n log n) | degree ≥ 1024 |

### Integer FFT via CRT

For exact integer polynomial multiplication, Tertius uses:
1. Three NTT-friendly primes: 998244353, 1224736769, 469762049
2. NTT multiplication in each field
3. Chinese Remainder Theorem reconstruction

This yields exact results for coefficients up to ~2^90.

### Sparse Multiplication

Geobucket algorithm with O(log n) amortized insertion time.

## Testing

```bash
# Run all tests (217 tests total)
cargo test

# Run property-based tests
cargo test proptests

# Run benchmarks
cargo bench
```

### Property-Based Testing

Using `proptest`, Tertius verifies:
- Ring axioms (associativity, commutativity, distributivity)
- Field axioms (multiplicative inverses)
- Polynomial evaluation homomorphism
- Karatsuba correctness vs schoolbook

## Comparing with SymPy

```bash
# Phase 1 NTT benchmarks
cd scripts
uv run python benchmark_sympy.py

# Phase 2 comprehensive benchmarks (multiplication, GCD, factorization, Gröbner)
uv run python benchmark_phase2.py
```

Sample output from Phase 2 benchmarks:
```
Polynomial Multiplication (NTT vs SymPy ZZ):
  Degree   Tertius NTT (ms)   SymPy ZZ (ms)   Speedup
      16              0.007           0.041      5.9x
      64              0.025           0.406     16.3x
     256              0.107           3.663     34.2x
    1024              0.512          26.510     51.8x

Note: Tertius NTT uses modular arithmetic (single prime).
      SymPy uses arbitrary precision integers.
```

## Crate Documentation

| Crate | Description |
|-------|-------------|
| `tertius-core` | Arena allocation, expression DAG, hash-consing |
| `tertius-integers` | `Integer`, `Rational`, `ModInt<P>` types |
| `tertius-rings` | `Ring`, `Field`, `EuclideanDomain` traits |
| `tertius-poly` | `DensePoly<R>`, `SparsePoly<R>`, FFT/NTT, sparse GCD |
| `tertius-simplify` | `egg`-based simplification, calculus rules |
| `tertius-linalg` | Sparse matrices (CSR), Block Wiedemann, Smith Normal Form |
| `tertius-factor` | Van Hoeij factorization, LLL, Berlekamp-Zassenhaus |
| `tertius-groebner` | M5GB Gröbner basis algorithm |
| `tertius` | Unified API facade |

## Contributing

Contributions are welcome! Areas of interest:

- [x] Gröbner basis computation (M5GB algorithm)
- [x] Factorization algorithms (Van Hoeij, LLL, Berlekamp-Zassenhaus)
- [x] Symbolic integration (differentiation/integration rules)
- [ ] Risch algorithm (full elementary function integration)
- [ ] Python bindings (PyO3)
- [ ] WASM target
- [ ] Multivariate factorization (Lecerf's algorithm)

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Acknowledgments

- [dashu](https://crates.io/crates/dashu) for arbitrary precision arithmetic
- [egg](https://crates.io/crates/egg) for equality saturation
- The Rust community for excellent tooling
