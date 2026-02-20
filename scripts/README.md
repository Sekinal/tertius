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
