# Phase 1: Foundational Infrastructure for Symbolic Integration

## Overview

This phase creates the foundational components needed for the Risch algorithm:
- Rational function type with canonical representation
- Extended GCD with Bézout coefficients
- Partial fraction decomposition
- Hermite reduction

## 1. tertius-rational-func Crate

### 1.1 RationalFunction Type

```rust
pub struct RationalFunction<K: Field> {
    numerator: DensePoly<K>,   // numerator polynomial
    denominator: DensePoly<K>, // monic, coprime with numerator
}
```

**Invariants:**
- Denominator is always monic (leading coefficient = 1)
- Numerator and denominator are coprime (gcd = 1)
- Zero is represented as 0/1

**Files:**
- `crates/tertius-rational-func/Cargo.toml`
- `crates/tertius-rational-func/src/lib.rs`
- `crates/tertius-rational-func/src/rational_func.rs`

### 1.2 Arithmetic Operations

Full field arithmetic for rational functions:
- Addition: a/b + c/d = (ad + bc) / bd → normalize
- Multiplication: (a/b) * (c/d) = ac / bd → normalize
- Division: (a/b) / (c/d) = ad / bc → normalize
- Negation: -(a/b) = (-a) / b

**File:** `crates/tertius-rational-func/src/arithmetic.rs`

### 1.3 Partial Fraction Decomposition

Given f = P/Q where deg(P) < deg(Q) and Q = q₁^{e₁} * ... * qₖ^{eₖ}:

```
P/Q = Σᵢ Σⱼ aᵢⱼ / qᵢʲ
```

**Algorithm:**
1. Factor denominator (or use squarefree decomposition)
2. For each factor qᵢ^{eᵢ}, compute coefficients via Euclidean algorithm

**File:** `crates/tertius-rational-func/src/partial_fractions.rs`

### 1.4 Hermite Reduction

Reduces ∫(A/D) to g + ∫(C/S) where:
- g is a rational function (no integration needed)
- S is squarefree (only simple poles remain)

**Algorithm:**
1. Compute squarefree decomposition: D = D₁ * D₂² * D₃³ * ...
2. Iteratively reduce higher powers using:
   - A/Dᵢ^j = B'/Dᵢ^{j-1} + C/Dᵢ^{j-1} for j > 1
   - Where B, C from: extended_gcd(Dᵢ, Dᵢ') then solve A = -B' * Dᵢ/gcd + C * Dᵢ'/gcd

**File:** `crates/tertius-rational-func/src/hermite.rs`

## 2. Extended GCD for Polynomials

### Location
`crates/tertius-poly/src/algorithms/gcd.rs`

### Function
```rust
pub fn poly_extended_gcd<F: Field>(
    a: &DensePoly<F>,
    b: &DensePoly<F>,
) -> (DensePoly<F>, DensePoly<F>, DensePoly<F>)
// Returns (gcd, s, t) where gcd = s*a + t*b
```

**Algorithm:** Standard extended Euclidean algorithm for polynomials.

## 3. Squarefree Decomposition

Required for Hermite reduction and partial fractions.

### Algorithm (Yun's algorithm)
Given polynomial f:
1. Compute g = gcd(f, f')
2. Squarefree part: f / g

For full decomposition into f = f₁ * f₂² * f₃³ * ...:
1. a₀ = gcd(f, f'), b₀ = f/a₀, c₀ = f'/a₀
2. Iterate: dᵢ = cᵢ - bᵢ', aᵢ₊₁ = gcd(bᵢ, dᵢ), bᵢ₊₁ = bᵢ/aᵢ₊₁, cᵢ₊₁ = dᵢ/aᵢ₊₁
3. Each aᵢ is squarefree factor with multiplicity i

**File:** `crates/tertius-poly/src/algorithms/squarefree.rs`

## 4. Implementation Order

1. ✅ Create tertius-rational-func crate structure
2. ✅ Implement RationalFunction type with normalize()
3. ✅ Implement arithmetic (Add, Sub, Mul, Div)
4. ⬜ Add poly_extended_gcd to tertius-poly
5. ⬜ Implement squarefree decomposition
6. ⬜ Implement partial fraction decomposition
7. ⬜ Implement Hermite reduction
8. ⬜ Add tests and benchmarks

## 5. Dependencies

```
tertius-rational-func
         ↓
    tertius-poly ←── poly_extended_gcd, squarefree
         ↓
    tertius-rings ←── Field trait
         ↓
   tertius-integers
```

## 6. Testing Strategy

- Unit tests for each operation
- Property tests:
  - normalize(a/b) is canonical
  - (a/b) * (b/a) = 1 for non-zero a/b
  - partial_fractions round-trips
  - Hermite: derivative of g equals A/D - C/S
