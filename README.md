# Tertius

**A high-performance Computer Algebra System in Rust**

Tertius is a modular CAS designed for speed and correctness. It leverages Rust's type system, modern algorithms, and parallel computation to achieve performance that significantly exceeds interpreted CAS implementations.

[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.75+-orange.svg)](https://www.rust-lang.org/)

## Performance Overview

Tertius consistently outperforms SymPy across polynomial algebra operations. The speedups come from:

1. **Compiled native code** vs Python interpreter overhead
2. **Cache-efficient data structures** (bit-packed monomials, arena allocation)
3. **Modern algorithms** (M5GB, Van Hoeij with LLL, FGLM)
4. **Parallel execution** via Rayon

### Benchmark Summary

| Category | Operation | Tertius | SymPy | Speedup |
|----------|-----------|---------|-------|---------|
| **Multiplication** | degree 1024 NTT | 0.51 ms | 33.6 ms | **66x** |
| **Factorization** | (x+1)^5 | 3.7 µs | 851 µs | **230x** |
| **Gröbner Basis** | katsura-3 | 35 µs | 1.1 ms | **31x** |
| **FGLM (Solving)** | quadratic system | 2.5 µs | 229 µs | **92x** |

*Benchmarks run on Linux x86_64. SymPy 1.13+, Python 3.12. Tertius compiled with `--release`.*

---

## Features

### Core Capabilities

- **Arbitrary Precision Arithmetic**: Integers and rationals with no size limits
- **Polynomial Algebra**: Dense and sparse representations with adaptive algorithms
- **Modular Arithmetic**: Compile-time and runtime prime fields (GF(p))
- **Algebraic Number Fields**: Q(α) for algebraic extensions
- **Symbolic Simplification**: Equality saturation via the `egg` library
- **Symbolic Integration**: Risch algorithm for elementary function integration

### Advanced Algorithms

- **Polynomial Factorization**: Van Hoeij algorithm with LLL lattice reduction
- **Gröbner Bases**: M5GB algorithm (F5 signatures + M4GB caching)
- **FGLM**: Gröbner basis conversion from grevlex to lex ordering
- **Sparse Linear Algebra**: Block Wiedemann, Smith Normal Form
- **Sparse Interpolation**: Ben-Or/Tiwari algorithm
- **Sparse GCD**: Hu-Monagan parallel algorithm
- **Risch Algorithm**: Complete rational function integration with Hermite reduction and Rothstein-Trager
- **Special Functions**: Polylogarithms, elliptic integrals, hypergeometric functions

---

## Detailed Benchmarks

### Polynomial Multiplication

NTT-based multiplication over finite fields vs SymPy's arbitrary-precision integer multiplication:

| Degree | Tertius NTT | SymPy ZZ | Speedup |
|--------|-------------|----------|---------|
| 16 | 7 µs | 43 µs | **6x** |
| 64 | 25 µs | 390 µs | **16x** |
| 256 | 107 µs | 3.7 ms | **35x** |
| 1024 | 512 µs | 33.6 ms | **66x** |

*Tertius uses NTT over a single prime field. SymPy uses arbitrary-precision integers (ZZ).*

### Polynomial Factorization

Van Hoeij factorization with LLL lattice reduction:

| Polynomial | Tertius | SymPy | Speedup |
|------------|---------|-------|---------|
| x⁴ - 1 | 25 µs | 426 µs | **17x** |
| x⁶ - 1 | 61 µs | 452 µs | **7x** |
| (x+1)⁵ | 3.7 µs | 851 µs | **230x** |
| x⁴ + 4 (Sophie Germain) | 38 µs | 657 µs | **17x** |

*Factorization over Z using Hensel lifting and LLL-based factor recombination.*

### Gröbner Basis Computation

M5GB algorithm (F5 signatures + matrix reduction):

| System | Variables | Tertius | SymPy | Speedup |
|--------|-----------|---------|-------|---------|
| cyclic-3 | 3 | 39 µs | 500 µs | **13x** |
| katsura-3 | 3 | 35 µs | 1.1 ms | **31x** |

*Grevlex ordering over GF(101).*

### FGLM Algorithm (System Solving)

Gröbner basis conversion from grevlex to lex for triangular solving:

| System | Solutions | Tertius | SymPy | Speedup |
|--------|-----------|---------|-------|---------|
| linear 2-var | 1 | 1.6 µs | 195 µs | **122x** |
| quadratic 2-var | 2 | 2.5 µs | 229 µs | **92x** |
| cyclic-2 | 2 | 3.0 µs | 170 µs | **57x** |

*Tertius: M5GB + FGLM. SymPy: groebner(order=lex).*

### Sparse Linear Algebra

CSR sparse matrix-vector multiply:

| Size | Tertius SpMV | Notes |
|------|--------------|-------|
| 10×10 | 68 ns | Bidiagonal |
| 50×50 | 327 ns | Bidiagonal |
| 100×100 | 660 ns | Bidiagonal |

### Polynomial GCD

Euclidean algorithm over finite fields:

| Input | Tertius | SymPy | Speedup |
|-------|---------|-------|---------|
| degree 5, common factor | 100 ns | 35 µs | **350x** |

---

## Architecture

```
tertius/
├── crates/
│   ├── tertius-core/         # Arena, expressions, hash-consing
│   ├── tertius-integers/     # Arbitrary precision (dashu wrapper)
│   ├── tertius-rings/        # Ring/Field traits, Z, Q, Z_p
│   ├── tertius-poly/         # Polynomial arithmetic + sparse algorithms
│   ├── tertius-simplify/     # egg-based simplification + calculus
│   ├── tertius-linalg/       # Sparse linear algebra (CSR, Wiedemann)
│   ├── tertius-factor/       # Polynomial factorization (Van Hoeij, LLL)
│   ├── tertius-groebner/     # Gröbner bases (M5GB algorithm)
│   ├── tertius-solve/        # FGLM, triangular solving
│   ├── tertius-rational-func/# Rational functions P(x)/Q(x)
│   ├── tertius-diffalg/      # Differential algebra, transcendental towers
│   ├── tertius-integrate/    # Risch algorithm, symbolic integration
│   ├── tertius-special-func/ # Polylogarithms, elliptic, hypergeometric
│   ├── tertius-series/       # Power series, Taylor expansion
│   ├── tertius-limits/       # Limit computation via Gruntz
│   └── tertius/              # Facade crate
├── benches/                  # Criterion benchmarks
├── examples/                 # Usage examples
└── scripts/                  # SymPy comparison scripts
```

### Design Principles

| Principle | Implementation |
|-----------|---------------|
| **Zero-cost abstractions** | Trait-based generics, monomorphization |
| **Cache efficiency** | Bit-packed monomials, arena allocation |
| **Parallelism** | Rayon-based data parallelism throughout |
| **Algorithm selection** | Automatic switching (schoolbook → Karatsuba → FFT) |

---

## Quick Start

```rust
use tertius::prelude::*;

// Polynomial arithmetic over rationals
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

---

## Algorithm Details

### Polynomial Multiplication

| Algorithm | Complexity | Threshold |
|-----------|------------|-----------|
| Schoolbook | O(n²) | degree < 32 |
| Karatsuba | O(n^1.58) | 32 ≤ degree < 1024 |
| NTT/FFT | O(n log n) | degree ≥ 1024 |

### Integer FFT via CRT

For exact integer polynomial multiplication:
1. Three NTT-friendly primes: 998244353, 1224736769, 469762049
2. NTT multiplication in each field
3. Chinese Remainder Theorem reconstruction

This yields exact results for coefficients up to ~2^90.

### Gröbner Basis (M5GB)

The M5GB algorithm combines:
- **F5 signatures**: Avoid redundant S-polynomial reductions
- **M4GB caching**: Incremental matrix construction
- **Parallel reduction**: Rayon-based row reduction

### FGLM Algorithm

Converts Gröbner bases from grevlex (fast to compute) to lex (triangular form):
1. Enumerate monomials in lex order
2. Compute normal forms w.r.t. grevlex basis
3. Detect linear dependencies → new lex basis elements

Complexity: O(nD³) where D is the dimension of K[x]/I.

### Risch Algorithm (Symbolic Integration)

The Risch algorithm is a complete decision procedure for elementary function integration:

**Rational Function Integration:**
1. **Hermite Reduction**: Decompose ∫(A/D) = g + ∫(C/S) where S is squarefree
2. **Rothstein-Trager**: Compute logarithmic part using resultants
   - res_t(A - t·D', D) gives the algebraic extension containing coefficients

**Transcendental Extensions:**
- **Logarithmic** (θ = log u): θ' = u'/u, integrate via polynomial reduction
- **Exponential** (θ = exp u): θ' = u'·θ, integrate via ansatz matching

**Non-Elementary Detection:**
- Uses Liouville's theorem to prove when no elementary antiderivative exists
- Recognizes canonical non-elementary forms: erf, Ei, li, Si, etc.

### Special Functions

| Function | Definition | Use Case |
|----------|------------|----------|
| Li_n(x) | Polylogarithm | ∫(log^n x)/x |
| erf(x) | Error function | ∫exp(-x²) |
| F(φ,k) | Elliptic integral 1st kind | ∫1/√(1-k²sin²θ) |
| E(φ,k) | Elliptic integral 2nd kind | ∫√(1-k²sin²θ) |
| ₂F₁(a,b;c;z) | Hypergeometric | General series solutions |

---

## Testing

```bash
# Run all tests
cargo test

# Run property-based tests
cargo test proptests

# Run benchmarks
cargo bench

# Compare with SymPy
cd scripts && uv run python benchmark_phase2.py
```

### Property-Based Testing

Using `proptest`, Tertius verifies:
- Ring axioms (associativity, commutativity, distributivity)
- Field axioms (multiplicative inverses)
- Polynomial evaluation homomorphism
- Karatsuba correctness vs schoolbook

---

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
| `tertius-solve` | FGLM algorithm, triangular solving |
| `tertius-rational-func` | Rational functions, partial fractions, Hermite reduction |
| `tertius-diffalg` | Differential algebra, transcendental tower K(θ₁)(θ₂)... |
| `tertius-integrate` | Risch algorithm, Rothstein-Trager, symbolic integration |
| `tertius-special-func` | Polylogarithms, elliptic integrals, hypergeometric ₂F₁ |
| `tertius-series` | Power series, Taylor expansion, asymptotic analysis |
| `tertius-limits` | Limit computation via Gruntz algorithm |
| `tertius` | Unified API facade |

---

## Contributing

Contributions are welcome! Areas of interest:

- [x] Gröbner basis computation (M5GB algorithm)
- [x] FGLM algorithm (grevlex → lex conversion)
- [x] Factorization algorithms (Van Hoeij, LLL, Berlekamp-Zassenhaus)
- [x] Symbolic integration (differentiation/integration rules)
- [x] Risch algorithm (full elementary function integration)
- [x] Special functions (polylogarithms, elliptic integrals, hypergeometric)
- [x] Non-integrability proofs (Liouville's theorem)
- [ ] Python bindings (PyO3)
- [ ] WASM target
- [ ] Multivariate factorization (Lecerf's algorithm)
- [ ] Differential equations (symbolic ODE solver)

---

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Acknowledgments

- [dashu](https://crates.io/crates/dashu) for arbitrary precision arithmetic
- [egg](https://crates.io/crates/egg) for equality saturation
- The Rust community for excellent tooling
