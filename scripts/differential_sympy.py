#!/usr/bin/env python3
"""Differential tests against SymPy for Tertius algorithms.

Current scope:
- Univariate integer polynomial factorization (degree multiset comparison)
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List

from sympy import Poly, ZZ, symbols


X = symbols("x")


@dataclass
class Failure:
    kind: str
    coeffs: List[int]
    details: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-degree", type=int, default=2)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--coeff-min", type=int, default=-8)
    parser.add_argument("--coeff-max", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=0.8)
    parser.add_argument("--max-failures", type=int, default=20)
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("target/debug/univariate_factor_cli"),
        help="Path to compiled univariate factor CLI",
    )
    parser.add_argument(
        "--build-binary",
        action="store_true",
        help="Build univariate factor CLI before running tests",
    )
    return parser.parse_args()


def build_binary() -> None:
    subprocess.run(
        ["cargo", "build", "-p", "tertius-factor", "--bin", "univariate_factor_cli"],
        check=True,
    )


def random_poly_coeffs(rng: random.Random, min_degree: int, max_degree: int, coeff_min: int, coeff_max: int) -> List[int]:
    degree = rng.randint(min_degree, max_degree)
    coeffs = [rng.randint(coeff_min, coeff_max) for _ in range(degree)]

    non_zero = [c for c in range(coeff_min, coeff_max + 1) if c != 0]
    coeffs.append(rng.choice(non_zero))
    if all(c == 0 for c in coeffs[:-1]):
        coeffs[0] = 1
    return coeffs


def parse_tertius_factor_output(raw: str) -> List[List[int]]:
    factors: List[List[int]] = []
    if not raw.strip():
        return factors

    for part in raw.strip().split(";"):
        if not part:
            continue
        coeffs = [int(item) for item in part.split(",") if item]
        while len(coeffs) > 1 and coeffs[-1] == 0:
            coeffs.pop()
        factors.append(coeffs)
    return factors


def factor_degrees_sympy(coeffs: List[int]) -> List[int]:
    poly = Poly(sum(c * X**i for i, c in enumerate(coeffs)), X, domain=ZZ)
    _, factors = poly.factor_list()
    degrees: List[int] = []
    for f, multiplicity in factors:
        degrees.extend([f.degree()] * multiplicity)
    degrees.sort()
    return degrees


def factor_degrees_tertius(binary: Path, coeffs: List[int], timeout: float) -> List[int]:
    coeff_arg = ",".join(str(c) for c in coeffs)
    proc = subprocess.run(
        [str(binary), coeff_arg],
        capture_output=True,
        text=True,
        timeout=timeout,
        check=True,
    )
    factors = parse_tertius_factor_output(proc.stdout)
    degrees = [len(f) - 1 for f in factors if not (len(f) == 1 and f[0] == 1)]
    degrees.sort()
    return degrees


def main() -> int:
    args = parse_args()

    if args.build_binary:
        build_binary()

    if not args.binary.exists():
        print(f"binary not found: {args.binary}", file=sys.stderr)
        print("run with --build-binary or build manually:", file=sys.stderr)
        print("cargo build -p tertius-factor --bin univariate_factor_cli", file=sys.stderr)
        return 2

    rng = random.Random(args.seed)

    failures: List[Failure] = []
    num_pass = 0
    num_timeout = 0
    num_exec_err = 0

    # Seed fixed regressions first.
    fixed_cases = [
        [1, -2, -8],
        [-3, -7, 0, 4],
        [-1, 1, 4, -8, -5, 8, 8, 8],
        [-1, 0, 4, -2, -4, 2, 5, -6],
        [2, -6, -1, -2, 4, 4, -8],
    ]

    for coeffs in fixed_cases:
        try:
            td = factor_degrees_tertius(args.binary, coeffs, args.timeout)
            sd = factor_degrees_sympy(coeffs)
            if td != sd:
                failures.append(
                    Failure("mismatch", coeffs, f"tertius={td} sympy={sd}")
                )
            else:
                num_pass += 1
        except subprocess.TimeoutExpired:
            num_timeout += 1
            failures.append(Failure("timeout", coeffs, "execution timed out"))
        except subprocess.CalledProcessError as exc:
            num_exec_err += 1
            stderr = (exc.stderr or "").strip()
            failures.append(Failure("exec_error", coeffs, stderr))

    for _ in range(args.cases):
        coeffs = random_poly_coeffs(
            rng,
            args.min_degree,
            args.max_degree,
            args.coeff_min,
            args.coeff_max,
        )
        try:
            td = factor_degrees_tertius(args.binary, coeffs, args.timeout)
            sd = factor_degrees_sympy(coeffs)
            if td != sd:
                failures.append(
                    Failure("mismatch", coeffs, f"tertius={td} sympy={sd}")
                )
            else:
                num_pass += 1
        except subprocess.TimeoutExpired:
            num_timeout += 1
            failures.append(Failure("timeout", coeffs, "execution timed out"))
        except subprocess.CalledProcessError as exc:
            num_exec_err += 1
            stderr = (exc.stderr or "").strip()
            failures.append(Failure("exec_error", coeffs, stderr))

        if len(failures) >= args.max_failures:
            break

    print(
        f"completed={num_pass + len(failures)} passed={num_pass} "
        f"timeouts={num_timeout} exec_errors={num_exec_err} failures={len(failures)}"
    )

    if failures:
        print("sample failures:")
        for failure in failures[: args.max_failures]:
            print(f"- kind={failure.kind} coeffs={failure.coeffs} details={failure.details}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
