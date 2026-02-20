//! LaTeX input/output helpers.

use smallvec::smallvec;
use tertius_core::arena::ExprArena;
use tertius_core::expr::{functions, ExprNode};
use tertius_core::handle::ExprHandle;
use tertius_units::Unit;

/// Converts an expression handle to LaTeX.
#[must_use]
pub fn expr_to_latex(arena: &ExprArena, expr: ExprHandle) -> String {
    match arena.get(expr) {
        ExprNode::Integer(n) => n.to_string(),
        ExprNode::Rational(num, den) => format!("\\frac{{{num}}}{{{den}}}"),
        ExprNode::Symbol(id) => arena.symbol_name(*id).unwrap_or("?").to_string(),
        ExprNode::Add(args) => {
            let mut out = String::new();
            for (i, arg) in args.iter().enumerate() {
                let piece = expr_to_latex(arena, *arg);
                if i == 0 {
                    out.push_str(&piece);
                    continue;
                }
                if let Some(rest) = piece.strip_prefix('-') {
                    out.push_str(" - ");
                    out.push_str(rest);
                } else {
                    out.push_str(" + ");
                    out.push_str(&piece);
                }
            }
            out
        }
        ExprNode::Mul(args) => args
            .iter()
            .map(|arg| {
                let p = expr_to_latex(arena, *arg);
                if matches!(arena.get(*arg), ExprNode::Add(_)) {
                    format!("\\left({p}\\right)")
                } else {
                    p
                }
            })
            .collect::<Vec<_>>()
            .join(" \\cdot "),
        ExprNode::Pow { base, exp } => {
            format!(
                "{{{}}}^{{{}}}",
                expr_to_latex(arena, *base),
                expr_to_latex(arena, *exp)
            )
        }
        ExprNode::Neg(arg) => format!("-{}", expr_to_latex(arena, *arg)),
        ExprNode::Div { num, den } => {
            format!(
                "\\frac{{{}}}{{{}}}",
                expr_to_latex(arena, *num),
                expr_to_latex(arena, *den)
            )
        }
        ExprNode::Function { id, args } => {
            let name = match *id {
                functions::SIN => "\\sin",
                functions::COS => "\\cos",
                functions::TAN => "\\tan",
                functions::EXP => "\\exp",
                functions::LN => "\\ln",
                functions::LOG10 => "\\log_{10}",
                functions::SQRT => "\\sqrt",
                functions::ABS => "\\left|",
                _ => "\\operatorname{f}",
            };
            if *id == functions::SQRT {
                format!("\\sqrt{{{}}}", expr_to_latex(arena, args[0]))
            } else if *id == functions::ABS {
                format!("\\left|{}\\right|", expr_to_latex(arena, args[0]))
            } else {
                format!("{name}\\left({}\\right)", expr_to_latex(arena, args[0]))
            }
        }
    }
}

/// Parses a limited linear-expression subset from LaTeX-like input.
///
/// Supported forms include terms like `2x - 3y + 5`.
pub fn parse_linear_latex(arena: &mut ExprArena, latex: &str) -> Result<ExprHandle, String> {
    let mut normalized = latex
        .replace(' ', "")
        .replace("\\cdot", "")
        .replace("\\times", "")
        .replace("\\,", "");
    if normalized.is_empty() {
        return Err("empty expression".to_string());
    }
    if !normalized.starts_with('-') {
        normalized = format!("+{normalized}");
    }

    let mut terms = Vec::new();
    let mut start = 0usize;
    let bytes = normalized.as_bytes();
    for i in 1..bytes.len() {
        if bytes[i] == b'+' || bytes[i] == b'-' {
            let term = &normalized[start..i];
            if !term.is_empty() {
                terms.push(parse_term(arena, term)?);
            }
            start = i;
        }
    }
    let last = &normalized[start..];
    if !last.is_empty() {
        terms.push(parse_term(arena, last)?);
    }

    if terms.is_empty() {
        return Err("could not parse any terms".to_string());
    }
    if terms.len() == 1 {
        Ok(terms[0])
    } else {
        Ok(arena.add(smallvec::SmallVec::from_vec(terms)))
    }
}

/// Parses a quantity in the form `3.5\,\mathrm{m}`.
pub fn parse_quantity_latex(latex: &str) -> Result<(f64, Unit), String> {
    let s = latex.replace(' ', "");
    let marker = "\\mathrm{";
    let pos = s
        .find(marker)
        .ok_or_else(|| "expected \\mathrm{...} unit marker".to_string())?;
    let value_str = s[..pos].trim_end_matches("\\,");
    let value = value_str
        .parse::<f64>()
        .map_err(|_| "invalid numeric value".to_string())?;
    let rest = &s[(pos + marker.len())..];
    let end = rest
        .find('}')
        .ok_or_else(|| "missing closing brace for unit".to_string())?;
    let symbol = &rest[..end];
    let unit = match symbol {
        "m" => Unit::meter(),
        "s" => Unit::second(),
        "kg" => Unit::kilogram(),
        "cm" => Unit::centimeter(),
        "N" => Unit::newton(),
        _ => return Err(format!("unsupported unit symbol: {symbol}")),
    };
    Ok((value, unit))
}

/// Formats a quantity as `value\,\mathrm{unit}`.
#[must_use]
pub fn format_quantity_latex(value: f64, unit: &Unit) -> String {
    format!("{}\\,\\mathrm{{{}}}", trim_float(value), unit.symbol)
}

fn parse_term(arena: &mut ExprArena, term: &str) -> Result<ExprHandle, String> {
    if let Some((num, den)) = parse_simple_fraction(term) {
        return Ok(arena.intern(ExprNode::Rational(num, den)));
    }

    let sign = if term.starts_with('-') { -1 } else { 1 };
    let body = term.trim_start_matches('+').trim_start_matches('-');
    if body.is_empty() {
        return Err("invalid empty term".to_string());
    }

    let split = body
        .find(|c: char| c.is_ascii_alphabetic() || c == '\\')
        .unwrap_or(body.len());
    let (coeff_part, symbol_part) = body.split_at(split);

    if symbol_part.is_empty() {
        let mut value = body
            .parse::<i64>()
            .map_err(|_| format!("invalid numeric term: {term}"))?;
        value *= i64::from(sign);
        return Ok(arena.integer(value));
    }

    let coeff = if coeff_part.is_empty() {
        i64::from(sign)
    } else {
        let parsed = coeff_part
            .parse::<i64>()
            .map_err(|_| format!("invalid coefficient in term: {term}"))?;
        parsed * i64::from(sign)
    };

    let (symbol_name, exp) = parse_symbol_and_exponent(symbol_part)?;
    let sym = arena.symbol(&symbol_name);
    let powered = if let Some(power) = exp {
        let p = arena.integer(power);
        arena.pow(sym, p)
    } else {
        sym
    };

    if coeff == 1 {
        Ok(powered)
    } else if coeff == -1 {
        Ok(arena.neg(powered))
    } else {
        let c = arena.integer(coeff);
        Ok(arena.mul(smallvec![c, powered]))
    }
}

fn parse_symbol_and_exponent(symbol_part: &str) -> Result<(String, Option<i64>), String> {
    let mut symbol = symbol_part.to_string();
    if let Some(stripped) = symbol.strip_prefix('\\') {
        symbol = stripped.to_string();
    }
    let mut exp = None;
    if let Some(idx) = symbol.find('^') {
        let base = symbol[..idx].to_string();
        let exp_str = symbol[(idx + 1)..].trim_matches('{').trim_matches('}');
        let parsed = exp_str
            .parse::<i64>()
            .map_err(|_| "invalid exponent".to_string())?;
        symbol = base;
        exp = Some(parsed);
    }
    if symbol.is_empty() {
        return Err("missing symbol name".to_string());
    }
    Ok((symbol, exp))
}

fn parse_simple_fraction(term: &str) -> Option<(i64, u64)> {
    let sign = if term.starts_with('-') { -1 } else { 1 };
    let body = term.trim_start_matches('+').trim_start_matches('-');
    let marker = "\\frac{";
    if !body.starts_with(marker) {
        return None;
    }
    let rest = &body[marker.len()..];
    let num_end = rest.find('}')?;
    let num = rest[..num_end].parse::<i64>().ok()?;
    let after_num = rest.get((num_end + 1)..)?;
    let after_slash = after_num.strip_prefix('/')?;
    let after_open = after_slash.strip_prefix('{')?;
    let den_end = after_open.find('}')?;
    let den = after_open[..den_end].parse::<u64>().ok()?;
    Some((num * i64::from(sign), den))
}

fn trim_float(value: f64) -> String {
    let s = format!("{value}");
    if s.ends_with(".0") {
        s.trim_end_matches(".0").to_string()
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;

    fn eval(arena: &ExprArena, expr: ExprHandle, vars: &HashMap<&str, f64>) -> f64 {
        match arena.get(expr) {
            ExprNode::Integer(n) => *n as f64,
            ExprNode::Rational(num, den) => *num as f64 / *den as f64,
            ExprNode::Symbol(id) => {
                let name = arena.symbol_name(*id).unwrap();
                *vars.get(name).unwrap()
            }
            ExprNode::Add(args) => args.iter().map(|h| eval(arena, *h, vars)).sum(),
            ExprNode::Mul(args) => args.iter().map(|h| eval(arena, *h, vars)).product(),
            ExprNode::Pow { base, exp } => eval(arena, *base, vars).powf(eval(arena, *exp, vars)),
            ExprNode::Neg(arg) => -eval(arena, *arg, vars),
            ExprNode::Div { num, den } => eval(arena, *num, vars) / eval(arena, *den, vars),
            ExprNode::Function { id, args } => {
                let x = eval(arena, args[0], vars);
                match *id {
                    functions::SIN => x.sin(),
                    functions::COS => x.cos(),
                    functions::EXP => x.exp(),
                    _ => panic!("unsupported function"),
                }
            }
        }
    }

    #[test]
    fn test_expr_to_latex_for_polynomial() {
        let mut arena = ExprArena::new();
        let x = arena.symbol("x");
        let two = arena.integer(2);
        let x2 = arena.pow(x, two);
        let one = arena.integer(1);
        let expr = arena.add(smallvec![x2, one]);
        let latex = expr_to_latex(&arena, expr);
        assert!(latex.contains("{x}^{2}"));
        assert!(latex.contains("1"));
    }

    #[test]
    fn test_parse_linear_latex_evaluates_correctly() {
        let mut arena = ExprArena::new();
        let expr = parse_linear_latex(&mut arena, "2x - 3y + 5").unwrap();
        let vars = HashMap::from([("x", 2.0), ("y", 1.0)]);
        assert!((eval(&arena, expr, &vars) - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_parse_and_format_quantity_latex() {
        let (value, unit) = parse_quantity_latex("3.5\\,\\mathrm{cm}").unwrap();
        assert!((value - 3.5).abs() < 1e-12);
        assert_eq!(unit.symbol, "cm");
        let formatted = format_quantity_latex(value, &unit);
        assert_eq!(formatted, "3.5\\,\\mathrm{cm}");
    }
}
