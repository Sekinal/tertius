# Tertius Physics-Student Roadmap

## Goal

Build Tertius into the best CAS for undergraduate physics workflows by focusing on:

1. Correct answers with domain-safe algebra.
2. Fast symbolic + numeric workflows for common physics problems.
3. Readable, step-by-step derivations suitable for study and homework.

## 6-Month Success Criteria

1. `>= 90%` correctness on a curated undergraduate physics corpus.
2. `<= 200 ms` median solve time for routine symbolic tasks.
3. Step-by-step derivations for `>= 70%` of solved problems.

## Current Gaps

1. Assumptions/domain semantics are not yet first-class across all symbolic operations.
2. No dedicated ODE/PDE pipeline for core physics equation classes.
3. No units and dimensional analysis engine.
4. Special-function support needs stronger AST-native, rule-based semantics.
5. No dedicated physics UX layer (derivation traces + textbook-style output).

## Implementation Plan

## Phase 1 (Weeks 1-4): Correctness Foundation

1. Extend differential testing beyond factorization to simplification, differentiation, and integration.
2. Add metamorphic/property testing for algebraic identities, edge cases, and branch behavior.
3. Create a physics regression corpus (mechanics, E&M, thermo, QM) and run it in CI.
4. Gate merges on correctness score and regression failures.

Deliverables:

1. Differential harnesses in `scripts/`.
2. CI workflow gates for correctness and stability.
3. Initial corpus with at least 200 physics-aligned expressions/problems.

## Phase 2 (Weeks 3-8): Assumptions + Units Core

1. Add assumption context in `tertius-core` (`real`, `positive`, `nonzero`, integer domains).
2. Thread assumptions through `tertius-simplify`, `tertius-diff`, and `tertius-integrate`.
3. Add a new `tertius-units` crate for SI units, conversion, and dimension checking.
4. Emit diagnostics for dimensional inconsistencies and unsafe transformations.

Deliverables:

1. Assumption-aware simplification/integration decisions.
2. Basic unit algebra and conversion.
3. Dimension-check errors integrated into high-level APIs.

## Phase 3 (Weeks 6-12): Physics Calculus Stack

1. Add grad/div/curl/Laplacian in Cartesian, cylindrical, and spherical coordinates.
2. Add symbolic scaffolding for line/surface/volume integrals.
3. Implement ODE core set: separable, first-order linear, second-order constant-coefficient.
4. Add Laplace-based solving for standard initial-value problems.

Deliverables:

1. Coordinate-aware vector calculus APIs.
2. ODE solver coverage for common textbook forms.
3. Symbolic + numerical fallback with verification hooks.

## Phase 4 (Weeks 10-18): Special Functions + Solver Depth

1. Prioritize Bessel, Legendre, Hermite, Laguerre, spherical harmonics.
2. Add simplification, derivative, and integral transformation rules for these families.
3. Extend `tertius-solve` and `tertius-linalg` for matrix exponential/eigen workflows.
4. Add PDE starter methods: separation of variables for heat/wave/Laplace in canonical geometries.

Deliverables:

1. Physics-special-function backbone.
2. Better linear-system and modal-analysis symbolic support.
3. PDE baseline workflows for core undergraduate scenarios.

## Phase 5 (Weeks 16-24): Physics-Focused UX

1. Build step-by-step derivation traces with rule provenance.
2. Add LaTeX input/output and textbook-style formatting.
3. Add a high-level "physics mode" API (`solve + assumptions + units + derivation`).
4. Publish correctness and performance dashboards.

Deliverables:

1. Derivation trace format and renderer.
2. Student-facing physics workflow entry points.
3. Public benchmark/correctness reporting.

## Crate-by-Crate Ownership

1. `tertius-core`: assumptions model and expression semantics.
2. `tertius-simplify`: assumption-safe canonicalization and rewrite safety.
3. `tertius-diff`: vector calculus and coordinate systems.
4. `tertius-integrate`: robust dispatch and calculus workflows.
5. `tertius-special-func`: physics special-function families.
6. `tertius-solve`: ODE/PDE symbolic solvers.
7. `tertius-linalg`: eigensystems, matrix exponentials, and stability.
8. `scripts/` + CI: differential tests, corpus checks, and performance gates.
9. `tertius-units` (new): units and dimensional analysis.
10. `tertius-physics` (new facade): end-to-end physics problem API.

## First 10 Implementation Tickets

1. Define assumption model and propagation hooks in `tertius-core`.
2. Implement a basic units type system and algebra in `tertius-units`.
3. Create physics corpus format and seed with 200 canonical problems.
4. Add differential harness for simplification vs SymPy.
5. Add differential harness for differentiation vs SymPy.
6. Add integration verification harness via differentiation checks.
7. Implement cylindrical/spherical grad/div/curl operations.
8. Implement first-order linear ODE solver.
9. Implement Laplace transform + inverse subset for IVPs.
10. Implement derivation trace format: `RuleApplied`, `Before`, `After`, `Conditions`.

## Prioritization Rule

If there is a tradeoff between breadth and reliability, prefer reliability for common physics-student problems. The target is trust and speed in daily use, not maximum theoretical coverage in the first pass.
