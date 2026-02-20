#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import sympy as sp


@dataclass
class Case:
    mode: str
    expr: str
    var: str = "x"
    assume: str = ""
    label: str = ""


TOKEN_RE = re.compile(r"\(|\)|[^\s()]+")


def tokenize(expr: str) -> list[str]:
    return TOKEN_RE.findall(expr)


def parse_sexpr(expr: str) -> Any:
    tokens = tokenize(expr)
    idx = 0

    def parse_one() -> Any:
        nonlocal idx
        if idx >= len(tokens):
            raise ValueError("unexpected end of expression")

        tok = tokens[idx]
        idx += 1
        if tok == "(":
            if idx >= len(tokens):
                raise ValueError("missing operator")
            op = tokens[idx]
            idx += 1
            args: list[Any] = []
            while idx < len(tokens) and tokens[idx] != ")":
                args.append(parse_one())
            if idx >= len(tokens):
                raise ValueError("missing closing ')'")
            idx += 1  # consume ')'
            return (op, args)
        if tok == ")":
            raise ValueError("unexpected ')'")
        return tok

    root = parse_one()
    if idx != len(tokens):
        raise ValueError("trailing tokens")
    return root


def sexpr_to_sympy(node: Any) -> sp.Expr:
    if isinstance(node, str):
        try:
            return sp.Integer(int(node))
        except ValueError:
            # Physics-student workflows are overwhelmingly over reals.
            return sp.Symbol(node, real=True)

    op, args = node
    sargs = [sexpr_to_sympy(a) for a in args]

    if op == "+":
        return sp.Add(*sargs)
    if op == "*":
        return sp.Mul(*sargs)
    if op == "-":
        if len(sargs) == 1:
            return -sargs[0]
        if len(sargs) == 2:
            return sargs[0] - sargs[1]
        raise ValueError("'-' requires one or two arguments")
    if op == "/":
        if len(sargs) != 2:
            raise ValueError("'/' requires two arguments")
        return sargs[0] / sargs[1]
    if op == "^":
        if len(sargs) != 2:
            raise ValueError("'^' requires two arguments")
        return sp.Pow(sargs[0], sargs[1])
    if op == "neg":
        if len(sargs) != 1:
            raise ValueError("'neg' requires one argument")
        return -sargs[0]
    if op == "sin":
        return sp.sin(sargs[0])
    if op == "cos":
        return sp.cos(sargs[0])
    if op == "tan":
        return sp.tan(sargs[0])
    if op == "exp":
        return sp.exp(sargs[0])
    if op == "ln":
        return sp.log(sargs[0])
    if op == "sqrt":
        return sp.sqrt(sargs[0])
    if op == "abs":
        return sp.Abs(sargs[0])

    raise ValueError(f"unsupported operator in output: {op}")


def parse_assumptions(spec: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    if not spec:
        return out
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        sym, tag = item.split(":", 1)
        out.setdefault(sym.strip(), set()).add(tag.strip().lower())
    return out


def numeric_equivalent(
    lhs: sp.Expr,
    rhs: sp.Expr,
    assumptions: dict[str, set[str]],
    rng: random.Random,
    samples: int = 10,
) -> bool:
    symbols = sorted(lhs.free_symbols.union(rhs.free_symbols), key=lambda s: s.name)
    if not symbols:
        try:
            return abs(complex(sp.N(lhs)) - complex(sp.N(rhs))) < 1e-8
        except Exception:
            return False

    valid = 0
    for _ in range(samples * 3):
        subs = {}
        for sym in symbols:
            tags = assumptions.get(sym.name, set())
            if "positive" in tags or "pos" in tags:
                val = rng.uniform(0.2, 2.5)
            else:
                val = rng.uniform(-2.5, 2.5)
                if ("nonzero" in tags or "nz" in tags) and abs(val) < 0.2:
                    val = 0.2 if val >= 0 else -0.2
            subs[sym] = val

        try:
            lv = complex(sp.N(lhs.subs(subs)))
            rv = complex(sp.N(rhs.subs(subs)))
        except Exception:
            continue

        if not (abs(lv.real) < 1e12 and abs(lv.imag) < 1e12):
            continue
        if not (abs(rv.real) < 1e12 and abs(rv.imag) < 1e12):
            continue

        valid += 1
        if abs(lv - rv) > 1e-6:
            return False
        if valid >= samples:
            return True

    try:
        return sp.simplify(lhs - rhs) == 0
    except Exception:
        return False


def run_cli(
    binary: Path,
    mode: str,
    expr: str,
    var: str | None,
    assume: str,
    timeout: float,
) -> tuple[str, list[str]]:
    cmd = [str(binary), mode, "--expr", expr]
    if var:
        cmd.extend(["--var", var])
    if assume:
        cmd.extend(["--assume", assume])

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "cli failed")

    line = proc.stdout.strip().splitlines()
    if not line:
        raise RuntimeError("empty cli output")
    parts = line[-1].split("\t")
    return parts[0], parts[1:]


def load_corpus(path: Path) -> list[Case]:
    with path.open("r", encoding="utf-8") as f:
        raw = json.load(f)

    out: list[Case] = []
    for item in raw:
        out.append(
            Case(
                mode=item["mode"],
                expr=item["expr"],
                var=item.get("var", "x"),
                assume=item.get("assume", ""),
                label=item.get("label", ""),
            )
        )
    return out


def random_integrate_case(rng: random.Random) -> Case:
    choice = rng.choice(["poly", "sin", "cos", "exp", "inv"])
    if choice == "poly":
        n = rng.randint(0, 5)
        coeff = rng.randint(1, 5)
        if n == 0:
            expr = str(coeff)
        elif coeff == 1:
            expr = f"(^ x {n})"
        else:
            expr = f"(* {coeff} (^ x {n}))"
        return Case(mode="integrate", expr=expr, var="x", label=f"rnd_poly_{n}")
    if choice == "sin":
        return Case(mode="integrate", expr="(sin x)", var="x", label="rnd_sin")
    if choice == "cos":
        return Case(mode="integrate", expr="(cos x)", var="x", label="rnd_cos")
    if choice == "exp":
        return Case(mode="integrate", expr="(exp x)", var="x", label="rnd_exp")
    return Case(
        mode="integrate",
        expr="(/ 1 x)",
        var="x",
        assume="x:positive,x:real",
        label="rnd_inv",
    )


def random_diff_case(rng: random.Random) -> Case:
    n = rng.randint(1, 5)
    coeff = rng.randint(1, 6)
    base = f"(* {coeff} (^ x {n}))"
    extra = rng.choice(["(sin x)", "(cos x)", "(exp x)", "x"])
    expr = f"(+ {base} {extra})"
    return Case(mode="diff", expr=expr, var="x", label="rnd_diff")


def random_simplify_case(rng: random.Random) -> Case:
    base = rng.choice(["x", "(sin x)", "(cos x)", "(exp x)", "(+ x 1)", "(* x x)"])
    template = rng.choice(
        [
            "(+ {b} 0)",
            "(+ 0 {b})",
            "(* {b} 1)",
            "(* 1 {b})",
            "(ln (exp x))",
            "(exp (ln x))",
        ]
    )
    expr = template.format(b=base)
    assume = "x:real" if expr == "(ln (exp x))" else ""
    return Case(mode="simplify", expr=expr, assume=assume, label="rnd_simplify")


def verify_case(binary: Path, case: Case, timeout: float, rng: random.Random) -> tuple[bool, str]:
    assumptions = parse_assumptions(case.assume)
    tag, payload = run_cli(binary, case.mode, case.expr, case.var, case.assume, timeout)

    if case.mode in {"simplify", "diff"}:
        if tag != "OK" or not payload:
            return False, f"{case.label} expected OK, got {tag} {payload}"
        out_expr = payload[0]
        lhs = sexpr_to_sympy(parse_sexpr(case.expr))
        rhs = sexpr_to_sympy(parse_sexpr(out_expr))
        if case.mode == "diff":
            sym = sp.Symbol(case.var, real=True)
            lhs = sp.diff(lhs, sym)
        if not numeric_equivalent(lhs, rhs, assumptions, rng):
            return False, f"{case.label} mismatch: in={case.expr} out={out_expr}"
        return True, ""

    if case.mode == "integrate":
        if tag != "OK" or not payload:
            return False, f"{case.label} expected OK, got {tag} {payload}"
        antideriv_expr = payload[0]
        sym_in = sexpr_to_sympy(parse_sexpr(case.expr))
        sym_antideriv = sexpr_to_sympy(parse_sexpr(antideriv_expr))
        sym = sp.Symbol(case.var, real=True)
        sym_diff = sp.diff(sym_antideriv, sym)
        if not numeric_equivalent(sym_in, sym_diff, assumptions, rng):
            return False, f"{case.label} mismatch: d/d{case.var}({antideriv_expr}) != {case.expr}"
        return True, ""

    return False, f"unknown mode: {case.mode}"


def build_binary() -> None:
    subprocess.run(
        ["cargo", "build", "--release", "--bin", "symbolic_cli"],
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Differential symbolic checks vs SymPy")
    parser.add_argument("--build-binary", action="store_true")
    parser.add_argument("--binary", default="target/release/symbolic_cli")
    parser.add_argument("--corpus", default="scripts/corpus/physics_undergrad_core.json")
    parser.add_argument("--cases", type=int, default=300)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--timeout", type=float, default=0.7)
    args = parser.parse_args()

    if args.build_binary:
        build_binary()

    binary = Path(args.binary)
    if not binary.exists():
        print(f"binary not found: {binary}", file=sys.stderr)
        return 2

    rng = random.Random(args.seed)
    corpus = load_corpus(Path(args.corpus))
    random_cases: list[Case] = []
    for _ in range(args.cases):
        random_cases.append(random_simplify_case(rng))
        random_cases.append(random_diff_case(rng))
        random_cases.append(random_integrate_case(rng))

    all_cases = corpus + random_cases

    completed = 0
    passed = 0
    failures = 0
    timeouts = 0
    exec_errors = 0
    fail_examples: list[str] = []

    for case in all_cases:
        completed += 1
        try:
            ok, detail = verify_case(binary, case, args.timeout, rng)
            if ok:
                passed += 1
            else:
                failures += 1
                if len(fail_examples) < 10:
                    fail_examples.append(detail)
        except subprocess.TimeoutExpired:
            timeouts += 1
            if len(fail_examples) < 10:
                fail_examples.append(f"{case.label} timeout")
        except Exception as exc:  # noqa: BLE001
            exec_errors += 1
            if len(fail_examples) < 10:
                fail_examples.append(f"{case.label} exec error: {exc}")

    print(
        "completed={} passed={} failures={} timeouts={} exec_errors={}".format(
            completed, passed, failures, timeouts, exec_errors
        )
    )
    for item in fail_examples:
        print(f"  - {item}")

    if failures or timeouts or exec_errors:
        return 1
    return 0


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parents[1])
    raise SystemExit(main())
