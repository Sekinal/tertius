use std::env;
use std::process;

use tertius_factor::univariate::van_hoeij_factor;
use tertius_integers::Integer;
use tertius_poly::dense::DensePoly;
use tertius_rings::integers::Z;

fn parse_coefficients(arg: &str) -> Result<Vec<Z>, String> {
    let mut coeffs = Vec::new();
    for part in arg.split(',') {
        let trimmed = part.trim();
        if trimmed.is_empty() {
            continue;
        }
        let value = Integer::from_str_radix(trimmed, 10)
            .map_err(|e| format!("invalid integer '{trimmed}': {e}"))?;
        coeffs.push(Z(value));
    }

    if coeffs.is_empty() {
        return Err("no coefficients provided".to_string());
    }

    Ok(coeffs)
}

fn main() {
    let Some(arg) = env::args().nth(1) else {
        eprintln!("usage: univariate_factor_cli \"c0,c1,...,cn\"");
        process::exit(2);
    };

    let coeffs = match parse_coefficients(&arg) {
        Ok(c) => c,
        Err(msg) => {
            eprintln!("{msg}");
            process::exit(2);
        }
    };

    let poly = DensePoly::new(coeffs);
    let result = van_hoeij_factor(&poly);

    let mut encoded_factors = Vec::new();
    for factor in result.factors {
        let encoded = factor
            .coeffs()
            .iter()
            .map(|z| z.0.to_string())
            .collect::<Vec<_>>()
            .join(",");
        encoded_factors.push(encoded);
    }

    println!("{}", encoded_factors.join(";"));
}
