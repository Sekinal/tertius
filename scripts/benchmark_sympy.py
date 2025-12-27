#!/usr/bin/env python3
"""Benchmark SymPy polynomial multiplication for comparison with Tertius."""

import time
import random
from typing import List, Tuple

from sympy import symbols, Poly, ZZ


def random_poly(degree: int, max_coeff: int = 100) -> List[int]:
    """Generate random polynomial coefficients."""
    return [random.randint(-max_coeff, max_coeff) for _ in range(degree + 1)]


def benchmark_sympy_poly_mul(degree: int, iterations: int = 10) -> Tuple[float, float]:
    """Benchmark SymPy polynomial multiplication."""
    x = symbols('x')

    # Generate random polynomials
    coeffs_a = random_poly(degree)
    coeffs_b = random_poly(degree)

    # Create SymPy polynomials
    p = Poly(coeffs_a[::-1], x, domain=ZZ)  # SymPy uses descending order
    q = Poly(coeffs_b[::-1], x, domain=ZZ)

    # Warmup
    _ = p * q

    # Time multiple iterations
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        result = p * q
        end = time.perf_counter()
        times.append(end - start)

    mean = sum(times) / len(times)
    std = (sum((t - mean) ** 2 for t in times) / len(times)) ** 0.5

    return mean, std


def main():
    print("SymPy Polynomial Multiplication Benchmark")
    print("=" * 50)
    print()

    sizes = [16, 64, 256, 1024]

    print(f"{'Degree':<10} {'Mean Time (ms)':<20} {'Std Dev (ms)':<15}")
    print("-" * 50)

    for size in sizes:
        mean, std = benchmark_sympy_poly_mul(size)
        print(f"{size:<10} {mean * 1000:>15.3f} ms   {std * 1000:>10.3f} ms")

    print()
    print("Comparison with Tertius NTT (modular arithmetic):")
    print("-" * 50)

    # NTT times from Tertius benchmarks (single prime, very fast)
    tertius_ntt_times = {
        16: 0.007,    # ~7 µs
        64: 0.025,    # ~25 µs
        256: 0.107,   # ~107 µs
        1024: 0.512,  # ~512 µs
    }

    print(f"{'Degree':<10} {'Tertius NTT (ms)':<18} {'SymPy ZZ (ms)':<15} {'Speedup':<10}")
    print("-" * 55)

    for size in sizes:
        mean, _ = benchmark_sympy_poly_mul(size, iterations=5)
        tertius_ms = tertius_ntt_times.get(size, 0)
        sympy_ms = mean * 1000
        if tertius_ms > 0:
            speedup = sympy_ms / tertius_ms
            print(f"{size:<10} {tertius_ms:>14.3f}      {sympy_ms:>12.3f}   {speedup:>8.1f}x")
        else:
            print(f"{size:<10} {'N/A':>14}      {sympy_ms:>12.3f}   {'N/A':>8}")

    print()
    print("Note: Tertius NTT works in modular arithmetic (single prime).")
    print("      SymPy uses arbitrary precision integers over ZZ.")
    print("      Speedup > 1 means Tertius is faster.")


if __name__ == "__main__":
    random.seed(42)  # For reproducibility
    main()
