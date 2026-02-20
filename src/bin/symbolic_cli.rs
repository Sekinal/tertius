use std::env;

use smallvec::smallvec;
use tertius_core::assumptions::{Assumption, AssumptionSet};
use tertius_core::{ExprArena, ExprHandle, ExprNode};
use tertius_diff::diff_with_assumptions;
use tertius_integrate::{integrate_with_options, IntegrationOptions, IntegrationResult};
use tertius_simplify::Simplifier;

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        return Err(usage());
    }

    let command = &args[1];
    let mut expr = None::<String>;
    let mut var = None::<String>;
    let mut assume = None::<String>;

    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--expr" => {
                i += 1;
                expr = args.get(i).cloned();
            }
            "--var" => {
                i += 1;
                var = args.get(i).cloned();
            }
            "--assume" => {
                i += 1;
                assume = args.get(i).cloned();
            }
            other => return Err(format!("unknown option: {other}\n{}", usage())),
        }
        i += 1;
    }

    let assumptions = parse_assumptions(assume.as_deref())?;

    match command.as_str() {
        "simplify" => {
            let input = expr.ok_or_else(usage)?;
            let simplifier = Simplifier::new();
            let out = simplifier.simplify_str_with_assumptions(&input, &assumptions)?;
            println!("OK\t{out}");
        }
        "diff" => {
            let input = expr.ok_or_else(usage)?;
            let var_name = var.ok_or_else(usage)?;
            let mut arena = ExprArena::new();
            let parsed = parse_sexpr_to_expr(&mut arena, &input)?;
            let var_handle = arena.symbol(&var_name);
            let out = diff_with_assumptions(&mut arena, parsed, var_handle, &assumptions);
            println!("OK\t{}", expr_to_sexpr(&arena, out));
        }
        "integrate" => {
            let input = expr.ok_or_else(usage)?;
            let var_name = var.ok_or_else(usage)?;
            let mut arena = ExprArena::new();
            let parsed = parse_sexpr_to_expr(&mut arena, &input)?;
            let var_handle = arena.symbol(&var_name);
            let options = IntegrationOptions {
                assumptions: Some(assumptions),
                ..IntegrationOptions::default()
            };
            let result = integrate_with_options(&mut arena, parsed, var_handle, options);
            emit_integral_result(&arena, result);
        }
        _ => return Err(format!("unknown command: {command}\n{}", usage())),
    }

    Ok(())
}

fn usage() -> String {
    "Usage:\n  symbolic_cli simplify --expr \"<sexpr>\" [--assume \"x:real,x:positive\"]\n  symbolic_cli diff --expr \"<sexpr>\" --var x [--assume \"x:real\"]\n  symbolic_cli integrate --expr \"<sexpr>\" --var x [--assume \"x:real,x:positive\"]".to_string()
}

fn parse_assumptions(spec: Option<&str>) -> Result<AssumptionSet, String> {
    let mut out = AssumptionSet::new();
    let Some(spec) = spec else {
        return Ok(out);
    };
    if spec.trim().is_empty() {
        return Ok(out);
    }

    for pair in spec.split(',') {
        let mut it = pair.split(':');
        let symbol = it
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| format!("invalid assumption item: {pair}"))?;
        let tag = it
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| format!("invalid assumption item: {pair}"))?;
        if it.next().is_some() {
            return Err(format!("invalid assumption item: {pair}"));
        }

        let assumption = match tag.to_ascii_lowercase().as_str() {
            "real" => Assumption::Real,
            "positive" | "pos" => Assumption::Positive,
            "nonzero" | "non_zero" | "nz" => Assumption::NonZero,
            "integer" | "int" => Assumption::Integer,
            _ => return Err(format!("unsupported assumption tag: {tag}")),
        };
        out.assume(symbol, assumption);
    }

    Ok(out)
}

fn emit_integral_result(arena: &ExprArena, result: IntegrationResult) {
    match result {
        IntegrationResult::Symbolic(s) => {
            println!("OK\t{}\t{}", expr_to_sexpr(arena, s.result), s.display);
        }
        IntegrationResult::SpecialFunction(s) => println!("SPECIAL\t{}", s.display),
        IntegrationResult::Elliptic(s) => println!("ELLIPTIC\t{}", s.display),
        IntegrationResult::Numerical(n) => println!("NUMERIC\t{}", n),
        IntegrationResult::NonElementary(ne) => println!("NON_ELEM\t{}", ne.reason),
        IntegrationResult::Unknown(u) => println!("UNKNOWN\t{}", u),
    }
}

fn parse_sexpr_to_expr(arena: &mut ExprArena, input: &str) -> Result<ExprHandle, String> {
    let tokens = tokenize(input);
    if tokens.is_empty() {
        return Err("empty expression".to_string());
    }
    let mut idx = 0usize;
    let expr = parse_expr_tokens(arena, &tokens, &mut idx)?;
    if idx != tokens.len() {
        return Err("trailing tokens in expression".to_string());
    }
    Ok(expr)
}

fn tokenize(input: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    for ch in input.chars() {
        match ch {
            '(' | ')' => {
                if !current.is_empty() {
                    out.push(current.clone());
                    current.clear();
                }
                out.push(ch.to_string());
            }
            c if c.is_whitespace() => {
                if !current.is_empty() {
                    out.push(current.clone());
                    current.clear();
                }
            }
            _ => current.push(ch),
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    out
}

fn parse_expr_tokens(
    arena: &mut ExprArena,
    tokens: &[String],
    idx: &mut usize,
) -> Result<ExprHandle, String> {
    let token = tokens
        .get(*idx)
        .ok_or_else(|| "unexpected end of input".to_string())?;
    if token == "(" {
        *idx += 1;
        let op = tokens
            .get(*idx)
            .ok_or_else(|| "missing operator".to_string())?
            .clone();
        *idx += 1;

        let mut args = Vec::<ExprHandle>::new();
        loop {
            let next = tokens
                .get(*idx)
                .ok_or_else(|| "missing closing ')'".to_string())?;
            if next == ")" {
                *idx += 1;
                break;
            }
            args.push(parse_expr_tokens(arena, tokens, idx)?);
        }
        build_op(arena, &op, &args)
    } else if token == ")" {
        Err("unexpected ')'".to_string())
    } else {
        *idx += 1;
        if let Ok(n) = token.parse::<i64>() {
            Ok(arena.integer(n))
        } else {
            Ok(arena.symbol(token))
        }
    }
}

fn build_op(arena: &mut ExprArena, op: &str, args: &[ExprHandle]) -> Result<ExprHandle, String> {
    match op {
        "+" => {
            if args.is_empty() {
                return Err("'+' requires at least one argument".to_string());
            }
            Ok(arena.add(args.to_vec()))
        }
        "*" => {
            if args.is_empty() {
                return Err("'*' requires at least one argument".to_string());
            }
            Ok(arena.mul(args.to_vec()))
        }
        "-" => match args {
            [a] => Ok(arena.neg(*a)),
            [a, b] => {
                let neg_b = arena.neg(*b);
                Ok(arena.add(smallvec![*a, neg_b]))
            }
            _ => Err("'-' supports one or two arguments".to_string()),
        },
        "/" => match args {
            [num, den] => Ok(arena.intern(ExprNode::Div {
                num: *num,
                den: *den,
            })),
            _ => Err("'/' requires exactly two arguments".to_string()),
        },
        "^" => match args {
            [base, exp] => Ok(arena.pow(*base, *exp)),
            _ => Err("'^' requires exactly two arguments".to_string()),
        },
        "neg" => match args {
            [a] => Ok(arena.neg(*a)),
            _ => Err("'neg' requires exactly one argument".to_string()),
        },
        "sin" => make_function(arena, tertius_core::expr::functions::SIN, args),
        "cos" => make_function(arena, tertius_core::expr::functions::COS, args),
        "tan" => make_function(arena, tertius_core::expr::functions::TAN, args),
        "exp" => make_function(arena, tertius_core::expr::functions::EXP, args),
        "ln" => make_function(arena, tertius_core::expr::functions::LN, args),
        "log10" => make_function(arena, tertius_core::expr::functions::LOG10, args),
        "sqrt" => make_function(arena, tertius_core::expr::functions::SQRT, args),
        "abs" => make_function(arena, tertius_core::expr::functions::ABS, args),
        _ => Err(format!("unsupported operator/function: {op}")),
    }
}

fn make_function(
    arena: &mut ExprArena,
    id: u32,
    args: &[ExprHandle],
) -> Result<ExprHandle, String> {
    if args.len() != 1 {
        return Err("functions currently require one argument".to_string());
    }
    Ok(arena.intern(ExprNode::Function {
        id,
        args: smallvec![args[0]],
    }))
}

fn expr_to_sexpr(arena: &ExprArena, expr: ExprHandle) -> String {
    match arena.get(expr) {
        ExprNode::Integer(n) => n.to_string(),
        ExprNode::Rational(num, den) => format!("(/ {num} {den})"),
        ExprNode::Symbol(id) => arena.symbol_name(*id).unwrap_or("?").to_string(),
        ExprNode::Add(args) => {
            let inner: Vec<String> = args.iter().map(|h| expr_to_sexpr(arena, *h)).collect();
            format!("(+ {})", inner.join(" "))
        }
        ExprNode::Mul(args) => {
            let inner: Vec<String> = args.iter().map(|h| expr_to_sexpr(arena, *h)).collect();
            format!("(* {})", inner.join(" "))
        }
        ExprNode::Pow { base, exp } => {
            format!(
                "(^ {} {})",
                expr_to_sexpr(arena, *base),
                expr_to_sexpr(arena, *exp)
            )
        }
        ExprNode::Neg(arg) => format!("(neg {})", expr_to_sexpr(arena, *arg)),
        ExprNode::Div { num, den } => {
            format!(
                "(/ {} {})",
                expr_to_sexpr(arena, *num),
                expr_to_sexpr(arena, *den)
            )
        }
        ExprNode::Function { id, args } => {
            let name = match *id {
                tertius_core::expr::functions::SIN => "sin",
                tertius_core::expr::functions::COS => "cos",
                tertius_core::expr::functions::TAN => "tan",
                tertius_core::expr::functions::EXP => "exp",
                tertius_core::expr::functions::LN => "ln",
                tertius_core::expr::functions::LOG10 => "log10",
                tertius_core::expr::functions::SQRT => "sqrt",
                tertius_core::expr::functions::ABS => "abs",
                _ => "f",
            };
            let inner: Vec<String> = args.iter().map(|h| expr_to_sexpr(arena, *h)).collect();
            if name == "f" {
                format!("(f{} {})", id, inner.join(" "))
            } else {
                format!("({name} {})", inner.join(" "))
            }
        }
    }
}
