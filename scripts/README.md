## Differential Testing

Run randomized differential checks against SymPy for univariate factorization:

```bash
uv run --project scripts python scripts/differential_sympy.py --build-binary --cases 500
```

Useful options:

```bash
uv run --project scripts python scripts/differential_sympy.py \
  --build-binary \
  --cases 2000 \
  --seed 42 \
  --timeout 0.5
```

This command exits with a non-zero status on mismatches, timeouts, or execution errors,
so it can be used directly in CI jobs.

## Differential Testing (Simplify / Diff / Integrate)

Run symbolic differential checks against SymPy using the `symbolic_cli` binary:

```bash
uv run --project scripts python scripts/differential_symbolic_sympy.py --build-binary --cases 200
```

This runner checks:

1. Simplification equivalence.
2. Differentiation equivalence.
3. Integration correctness by differentiating returned antiderivatives.

It uses a curated physics-oriented corpus at `scripts/corpus/physics_undergrad_core.json`
and randomized cases. The script exits non-zero on failures/timeouts/execution errors.
